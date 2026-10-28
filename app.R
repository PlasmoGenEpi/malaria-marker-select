library(shiny)
library(DT)
library(readr)
library(dplyr)

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

# Helper to load the default marker table bundled with the app
load_default_markers <- function() {
  path <- file.path("data", "principal_resistance_marker_info_table.tsv")
  read_tsv(path, show_col_types = FALSE, progress = FALSE)
}

load_who_summary <- function() {
  path <- file.path("data", "who_compendium_marker_summary.tsv")
  read_tsv(path, show_col_types = FALSE, progress = FALSE) %>%
    mutate(
      gene = as.character(.data$gene),
      aa_position = as.numeric(.data$aa_position),
      in_who_compendium = TRUE
    )
}

load_who_long <- function() {
  path <- file.path("data", "who_compendium_annotations_long.tsv")
  read_tsv(path, show_col_types = FALSE, progress = FALSE) %>%
    mutate(
      gene = as.character(.data$gene),
      aa_position = as.numeric(.data$aa_position),
      drug = as.character(.data$drug),
      class_bucket = as.character(.data$class_bucket)
    )
}

who_summary <- load_who_summary()
who_long <- load_who_long()
who_drugs <- sort(unique(who_long$drug))
who_class_buckets <- c("Validated", "Candidate", "Potential")

# Read a BED / BED-like TSV into chrom, start, end (0-based half-open assumed).
read_bed_intervals <- function(path) {
  raw_lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  raw_lines <- raw_lines[!grepl("^(track|browser)\\b", raw_lines, ignore.case = TRUE)]
  raw_lines <- raw_lines[nzchar(trimws(raw_lines))]
  if (length(raw_lines) == 0) {
    stop("The BED file is empty.")
  }

  first <- raw_lines[[1]]
  has_header <- grepl("^#?chrom(\\t|,|$)", first, ignore.case = TRUE) ||
    grepl("^#chrom", first, ignore.case = TRUE)
  body <- if (has_header) raw_lines[-1] else raw_lines
  if (length(body) == 0) {
    stop("The BED file has a header but no intervals.")
  }

  tmp <- tempfile(fileext = ".bed")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(body, tmp)

  bed <- read_tsv(
    tmp,
    col_names = FALSE,
    col_types = cols(.default = col_character()),
    progress = FALSE,
    show_col_types = FALSE
  )
  if (ncol(bed) < 3) {
    stop("BED files need at least three columns: chrom, start, end.")
  }

  tibble(
    chrom = as.character(bed[[1]]),
    start = suppressWarnings(as.integer(bed[[2]])),
    end = suppressWarnings(as.integer(bed[[3]]))
  ) %>%
    filter(
      !is.na(.data$chrom), .data$chrom != "",
      !is.na(.data$start), !is.na(.data$end),
      .data$end > .data$start
    )
}

# Markers whose [start, stop) overlaps any panel interval on the same chrom.
markers_overlapping_bed <- function(markers_df, bed_df) {
  if (is.null(bed_df) || nrow(bed_df) == 0) {
    return(markers_df)
  }
  m <- markers_df %>%
    mutate(
      chrom = as.character(.data$chrom),
      start = as.integer(.data$start),
      stop = as.integer(.data$stop)
    )
  keep <- vapply(seq_len(nrow(m)), function(i) {
    chrom_i <- m$chrom[[i]]
    start_i <- m$start[[i]]
    stop_i <- m$stop[[i]]
    same <- bed_df$chrom == chrom_i
    if (!any(same)) {
      return(FALSE)
    }
    b <- bed_df[same, , drop = FALSE]
    any(start_i < b$end & stop_i > b$start)
  }, logical(1))
  markers_df[keep, , drop = FALSE]
}

# STAVE-style group id without amino acids: GENE:pos_pos;GENE:pos
stave_group_name <- function(df) {
  if (is.null(df) || nrow(df) == 0) {
    return("")
  }
  parts <- df %>%
    mutate(
      gene_id = as.character(.data$gene_id),
      aa_position = as.integer(.data$aa_position)
    ) %>%
    distinct(.data$gene_id, .data$aa_position) %>%
    arrange(.data$gene_id, .data$aa_position) %>%
    group_by(.data$gene_id) %>%
    summarise(
      positions = paste(.data$aa_position, collapse = "_"),
      .groups = "drop"
    ) %>%
    mutate(label = paste0(.data$gene_id, ":", .data$positions))
  paste(parts$label, collapse = ";")
}

annotate_with_who <- function(df) {
  df %>%
    mutate(
      gene = as.character(.data$gene),
      aa_position = as.numeric(.data$aa_position)
    ) %>%
    select(-any_of(c(
      "in_who_compendium", "who_drugs", "who_class_buckets",
      "who_classifications", "who_alterations"
    ))) %>%
    left_join(who_summary, by = c("gene", "aa_position")) %>%
    mutate(
      in_who_compendium = ifelse(is.na(.data$in_who_compendium), FALSE, TRUE),
      who_drugs = ifelse(is.na(.data$who_drugs), "", .data$who_drugs),
      who_class_buckets = ifelse(is.na(.data$who_class_buckets), "", .data$who_class_buckets),
      who_classifications = ifelse(is.na(.data$who_classifications), "", .data$who_classifications),
      who_alterations = ifelse(is.na(.data$who_alterations), "", .data$who_alterations)
    )
}

ui <- fluidPage(
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "app.css"),
    tags$meta(name = "viewport", content = "width=device-width, initial-scale=1")
  ),
  div(
    class = "app-header",
    tags$img(src = "logo.png", alt = "PlasmoGenEpi logo"),
    div(
      class = "app-header-text",
      h1("Malaria Marker Selector"),
      p(class = "app-tagline", "Filter, select, and download markers genomic locations")
    )
  ),
  sidebarLayout(
    sidebarPanel(
      helpText(
        "Markers include those from the WHO compendium, plus a curated selection of additional markers of interest."
      ),
      tags$details(
        class = "sidebar-details sidebar-details-filters",
        open = NA,
        tags$summary("Filters"),
        h4("WHO Compendium"),
        helpText(
          "Limit markers using WHO compendium annotations."
        ),
        checkboxInput(
          "who_only",
          "WHO compendium only",
          value = FALSE
        ),
        selectInput(
          "who_drug",
          "Drug(s)",
          choices = who_drugs,
          selected = character(0),
          multiple = TRUE
        ),
        selectInput(
          "who_class",
          "Classification",
          choices = who_class_buckets,
          selected = character(0),
          multiple = TRUE
        ),
        actionButton("clear_who_filters", "Clear WHO filters"),
        tags$hr(),
        h4("Panel BED"),
        helpText(
          "Upload a BED file (chrom, start, end) to only keep markers that are covered by your panel."
        ),
        fileInput(
          "panel_bed_file",
          "Upload panel BED",
          accept = c(".bed", ".tsv", ".bed.tsv", "text/tab-separated-values", "text/plain")
        ),
        conditionalPanel(
          condition = "output.panel_bed_loaded",
          actionButton("clear_panel_bed", "Clear panel BED")
        ),
        uiOutput("panel_bed_status")
      ),
      uiOutput("who_filter_status"),
      tags$details(
        class = "sidebar-details sidebar-details-columns",
        tags$summary("Columns"),
        helpText("Choose extra columns to show in the table."),
        uiOutput("extra_cols_ui")
      ),
      tags$details(
        class = "sidebar-details sidebar-details-custom",
        tags$summary("Custom Markers"),
        helpText("Replace the built-in marker table with your own TSV."),
        fileInput(
          "marker_file",
          "Upload marker TSV",
          accept = c(".tsv", "text/tab-separated-values")
        ),
        actionButton("reset_data", "Use built-in table")
      ),
      width = 3
    ),
    mainPanel(
      width = 9,
      div(
        class = "main-panel-card",
        tabsetPanel(
        tabPanel(
          "Markers",
          h3("Markers"),
          helpText(
            "Select individual markers to download their genomic coordinates. Use the sidebar to filter first if needed."
          ),
          div(
            style = "margin-bottom: 10px; display: flex; gap: 8px; flex-wrap: wrap;",
            actionButton("select_all_markers", "Select all"),
            actionButton("clear_marker_selection", "Clear selection")
          ),
          DTOutput("marker_table"),
          tags$hr(),
          h3("Download"),
          helpText(
            "Downloads a BED file with ", code("#chrom"), ", ", code("start"), ", ", code("end"), ", ",
            code("gene_id"), ", ", code("aa_position"), ", ", code("strand"), ", plus any optional columns you pick."
          ),
          uiOutput("download_cols_ui"),
          checkboxInput("show_markers_preview", "Show preview", value = FALSE),
          conditionalPanel(
            condition = "input.show_markers_preview == true",
            h4("Preview"),
            DTOutput("markers_download_preview")
          ),
          downloadButton("download_markers", "Download selected markers"),
          helpText(
            "Next: open ", strong("Marker groups"), " to assign markers to named groups if required."
          )
        ),
        tabPanel(
          "Marker groups",
          h3("Marker groups"),
          helpText(
            "Create marker groups from markers you selected on the Markers tab."
          ),
          uiOutput("grouping_placeholder_ui"),
          DTOutput("grouping_marker_table"),
          tags$hr(),
          h3("Create a group"),
          helpText(
            "Select rows above, name the group or check Autogenerate name, then add.",
            " A generated name is the gene id and positions, separated like a STAVE id, without amino acids.",
            " Markers can be in more than one group."
          ),
          div(
            class = "group-name-row",
            textInput(
              "group_name",
              "Group name",
              value = "",
              placeholder = "e.g. PF3D7_0417200:51_108"
            ),
            checkboxInput("autogenerate_group_name", "Autogenerate name", value = FALSE)
          ),
          div(
            style = "display: flex; gap: 8px; flex-wrap: wrap; margin-bottom: 12px;",
            actionButton("add_to_group", "Add to group"),
            actionButton("clear_groups", "Clear all groups")
          ),
          h4("Groups preview"),
          uiOutput("groups_download_cols_ui"),
          DTOutput("group_table"),
          downloadButton("download_groups", "Download marker groups")
        ),
        tabPanel(
          "Instructions",
          h3("How to use this app"),
          p(
            "Filter and select drug-resistance markers, download their genomic coordinates,",
            " and group any you want to use together. A marker is one locus.",
            " A group is several loci used together."
          ),
          p(
            "This app was designed to create inputs for the ",
            tags$a(
              href = "https://nf-co.re/plasmodiumdrugres",
              target = "_blank",
              rel = "noopener noreferrer",
              .noWS = "after",
              "plasmodiumdrugres pipeline"
            ),
            "."
          ),
          tags$ol(
            class = "instruction-steps",
            tags$li(
              strong("Filter, if you want a shorter list."),
              " Use the sidebar to keep markers in the WHO Compendium, their drug association, or classification",
              " (Validated, Candidate, or Potential). Upload a BED file to only keep markers that are covered by your panel/ regions of interest.",
            ),
            tags$li(
              strong("Select markers."),
              " On ", strong("Markers", .noWS = "after"), ", search or sort, then select the rows you want."
            ),
            tags$li(
              strong("Download coordinates."),
              " The BED file has ",
              code("#chrom", .noWS = "after"), ", ",
              code("start", .noWS = "after"), ", ",
              code("end", .noWS = "after"), ", ",
              code("gene_id", .noWS = "after"), ", ",
              code("aa_position", .noWS = "after"), ", and ",
              code("strand", .noWS = "after"), ". ",
              code("gene"), " is included by default. Optionally, add extra columns including WHO Compendium annotations."
            ),
            tags$li(
              strong("Make groups for multi-locus sets."),
              " On ", strong("Marker groups", .noWS = "after"),
              ", start from the markers you selected on ", strong("Markers", .noWS = "after"),
              ". Select a subset, name the group, and add it, then repeat for another subset.",
              " For example, you could put ", code("crt", .noWS = "after"),
              " mutations you selected on ", strong("Markers", .noWS = "after"), " in one group and ", code("dhfr"), " with ", code("dhps"),
              " in another. The download has ",
              code("group_id", .noWS = "after"), ", ",
              code("gene_id", .noWS = "after"), ", and ",
              code("aa_position", .noWS = "after"), "."
            )
          ),
          h3("Custom markers"),
          p(
            "To use your own table, upload a TSV with ",
            code("gene_id", .noWS = "after"), ", ",
            code("gene", .noWS = "after"), ", ",
            code("aa_position", .noWS = "after"), ", ",
            code("chrom", .noWS = "after"), ", ",
            code("start", .noWS = "after"), ", ",
            code("stop", .noWS = "after"), ", and ",
            code("strand", .noWS = "after"),
            ". WHO notes are added when the gene and amino-acid position match."
          ),
          p(
            "Column definitions are on the ", strong("Column guide"), " tab."
          )
        ),
        tabPanel(
          "Column guide",
          h3("Marker table columns"),
          p(
            "The Markers table is one row per amino‑acid site (with genomic coordinates).",
            " Use the sidebar ", strong("Columns"), " control to show extra fields.",
            " Core location columns come from the built‑in marker table; WHO columns are",
            " joined from the Compendium and ", em("collapsed across all WHO rows"),
            " that include that site."
          ),
          h4("Core columns"),
          tags$ul(
            tags$li(strong(code("gene_id")), " — PlasmoDB gene identifier used in downloads (e.g. ", code("PF3D7_0709000"), ") "),
            tags$li(strong(code("gene")), " — short gene name (", code("crt"), ", ", code("k13"), ", ", code("dhfr"), ", …)."),
            tags$li(strong(code("aa_position")), " — amino‑acid position in the gene."),
            tags$li(strong(code("chrom")), " / ", strong(code("start")), " / ", strong(code("stop")), " / ", strong(code("strand")), " — 0‑based half‑open genomic interval for that codon (BED‑like). Downloads rename ", code("stop"), " to ", code("end"), " and prefix chrom as ", code("#chrom"), ".")
          ),
          h4("WHO Compendium Annotation Columns"),
          p(
            "In the WHO Compendium Excel workbook, each row is usually one drug × one alteration (mutation or haplotype).",
            " Here the table is keyed by ", code("gene"), " + ", code("aa_position"), ", so every Excel row that",
            " mentions that site is merged into one marker row. Multi‑value fields are joined with ",
            code("; "), " (sorted, unique). See the ", strong("WHO Compendium field guide"), " for more details."
          ),
          tags$ul(
            tags$li(
              strong(code("in_who_compendium")),
              " — ", code("TRUE"), " if this site appears in at least one WHO AA alteration; otherwise ", code("FALSE"), "."
            ),
            tags$li(
              strong(code("who_drugs")),
              " — all drugs linked to this site across WHO rows (e.g. ", code("Chloroquine; Piperaquine"), ")."
            ),
            tags$li(
              strong(code("who_class_buckets")),
              " — simplified evidence classes for filtering: ", code("Validated"), ", ",
              code("Candidate"), ", and/or ", code("Potential"), ". A site can have more than one."
            ),
            tags$li(
              strong(code("who_classifications")),
              " — fuller WHO classification text from the Excel (e.g. \"Validated marker (Epi., Lab. & Clin.)\"),",
              " also collapsed with ", code("; "), " when several apply."
            ),
            tags$li(
              strong(code("who_alterations")),
              " — WHO ", em("Alteration(s)"), " strings that include this site. Unlike the Excel,",
              " this is ", strong("not"), " a single haplotype for the row: it lists every distinct",
              " alteration/haplotype string from the Compendium that contains this amino‑acid position,",
              " separated by ", code("; "), ". Within each string, WHO usually uses ", code("+"),
              " for co‑occurring sites in one gene and ", code("&"), " to join partners across genes",
              " (commonly ", code("dhfr"), " & ", code("dhps"), "). Example: a ", code("dhfr"), " 51 row may show",
              " both a single‑site change (", code("N51I"), ") and multi‑site haplotypes",
              " (", code("N51I + C59R + S108N & A437G + K540E"), "). Use these strings as a guide when building",
              " multi‑locus groups; they are not auto‑converted into Marker groups."
            )
          ),
          p(
            "Amplifications and other non‑AA WHO entries are not expanded to BED coordinates here."
          ),
          h3("WHO field guide"),
          p(
            "From the WHO Compendium user guide (Dec 2025, v1.0). Classifications reflect",
            " laboratory, clinical, and/or genetic epidemiological evidence."
          ),
          tags$ul(
            tags$li(strong("Gene(s):"), " gene(s) carrying the alteration."),
            tags$li(strong("Alteration(s):"), " mutation, haplotype, or amplification (shown in the app as ", code("who_alterations"), "; see above for collapsing)."),
            tags$li(strong("Classification:"), " Validated, Candidate, or Potential based on evidence types."),
            tags$li(strong("Evidence A1–C1:"), " 1 = threshold met — lab (A1 transfection; A2 isolate/strain assays; A3 genetic-cross progeny); clinical (B1 delayed clearance; B2 recrudescence; B3 post-treatment selection); epi (C1 prevalence)."),
            tags$li(strong("References (PMIDs):"), " in the Excel workbook only.")
          ),
          p(
            strong("Gene abbreviations:"),
            " ", code("Pfcrt"), " chloroquine resistance transporter; ",
            code("Pfcytb"), " cytochrome b; ",
            code("Pfdhfr"), " dihydrofolate reductase; ",
            code("Pfdhps"), " dihydropteroate synthase; ",
            code("Pfgch1"), " GTP-cyclohydrolase I; ",
            code("Pfk13"), " Kelch 13; ",
            code("Pfmdr1"), " multidrug resistance gene 1; ",
            code("Pfpm2/3"), " plasmepsin 2/3; ",
            code("Pfcoronin"), " coronin.",
            " The app table uses short names (", code("crt"), ", ", code("k13"), ", …)."
          )
        )
      )
      )
    )
  )
)

server <- function(input, output, session) {
  # Reactive marker table, starting from the built-in data
  markers <- reactiveVal(load_default_markers())

  observeEvent(input$marker_file, {
    req(input$marker_file)
    df <- read_tsv(input$marker_file$datapath, show_col_types = FALSE, progress = FALSE)
    required_cols <- c("gene_id", "gene", "aa_position", "chrom", "start", "stop", "strand")
    missing <- setdiff(required_cols, names(df))
    if (length(missing) > 0) {
      showModal(modalDialog(
        title = "Problem with uploaded file",
        paste(
          "The uploaded file is missing these required columns:",
          paste(missing, collapse = ", ")
        ),
        easyClose = TRUE
      ))
    } else {
      markers(df)
    }
  })

  observeEvent(input$reset_data, {
    markers(load_default_markers())
  })

  observeEvent(input$clear_who_filters, {
    updateCheckboxInput(session, "who_only", value = FALSE)
    updateSelectInput(session, "who_drug", selected = character(0))
    updateSelectInput(session, "who_class", selected = character(0))
  })

  panel_bed <- reactiveVal(NULL)
  output$panel_bed_loaded <- reactive({
    !is.null(panel_bed())
  })
  outputOptions(output, "panel_bed_loaded", suspendWhenHidden = FALSE)

  observeEvent(input$panel_bed_file, {
    req(input$panel_bed_file)
    tryCatch({
      bed <- read_bed_intervals(input$panel_bed_file$datapath)
      if (nrow(bed) == 0) {
        showModal(modalDialog(
          title = "Problem with panel BED",
          "No valid intervals were found. Check that the file has chrom, start, and end columns.",
          easyClose = TRUE
        ))
        return()
      }
      panel_bed(bed)
    }, error = function(e) {
      showModal(modalDialog(
        title = "Problem with panel BED",
        conditionMessage(e),
        easyClose = TRUE
      ))
    })
  })

  observeEvent(input$clear_panel_bed, {
    panel_bed(NULL)
  })

  annotated_markers <- reactive({
    annotate_with_who(markers())
  })

  filtered_markers <- reactive({
    df <- annotated_markers()
    drugs <- input$who_drug %||% character(0)
    classes <- input$who_class %||% character(0)
    who_only <- isTRUE(input$who_only)

    if (who_only || length(drugs) > 0 || length(classes) > 0) {
      key_df <- who_long
      if (length(drugs) > 0) {
        key_df <- key_df %>% filter(.data$drug %in% drugs)
      }
      if (length(classes) > 0) {
        key_df <- key_df %>% filter(.data$class_bucket %in% classes)
      }
      keep_keys <- key_df %>%
        distinct(.data$gene, .data$aa_position)
      df <- df %>% inner_join(keep_keys, by = c("gene", "aa_position"))
    }

    bed <- panel_bed()
    if (!is.null(bed)) {
      df <- markers_overlapping_bed(df, bed)
    }

    df
  })

  output$who_filter_status <- renderUI({
    n_all <- nrow(annotated_markers())
    n_show <- nrow(filtered_markers())
    div(
      class = "marker-count-status",
      strong(sprintf("Showing %d of %d markers", n_show, n_all))
    )
  })

  output$panel_bed_status <- renderUI({
    bed <- panel_bed()
    if (is.null(bed)) {
      return(helpText("No panel BED loaded."))
    }
    helpText(sprintf(
      "Panel BED loaded: %d intervals on %d chromosomes.",
      nrow(bed),
      dplyr::n_distinct(bed$chrom)
    ))
  })

  # Default and extra columns for display
  base_display_cols <- c("gene_id", "gene", "aa_position")
  who_display_cols <- c("in_who_compendium", "who_drugs", "who_class_buckets")

  output$extra_cols_ui <- renderUI({
    df <- annotated_markers()
    all_cols <- names(df)
    extra_choices <- setdiff(all_cols, base_display_cols)
    default_selected <- intersect(who_display_cols, extra_choices)
    checkboxGroupInput(
      "extra_cols",
      "Select additional columns to display",
      choices = extra_choices,
      selected = default_selected
    )
  })

  output$marker_table <- renderDT({
    df <- filtered_markers()
    display_cols <- unique(c(base_display_cols, input$extra_cols %||% character(0)))
    display_cols <- intersect(display_cols, names(df))
    datatable(
      df[, display_cols, drop = FALSE],
      selection = "multiple",
      filter = "top",
      width = "100%",
      options = list(
        pageLength = 25,
        autoWidth = FALSE,
        scrollX = TRUE,
        scrollCollapse = TRUE
      )
    )
  })

  marker_table_proxy <- dataTableProxy("marker_table")

  observeEvent(input$select_all_markers, {
    # Prefer rows that remain after DataTables search/column filters;
    # fall back to every row in the currently WHO-filtered table.
    rows <- input$marker_table_rows_all
    if (is.null(rows) || length(rows) == 0) {
      rows <- seq_len(nrow(filtered_markers()))
    }
    selectRows(marker_table_proxy, rows)
  })

  observeEvent(input$clear_marker_selection, {
    selectRows(marker_table_proxy, NULL)
  })

  # Extra columns for BED-like download
  output$download_cols_ui <- renderUI({
    df <- annotated_markers()
    base_bed_cols <- c("chrom", "start", "stop", "strand", "gene_id", "aa_position")
    extra_choices <- setdiff(names(df), base_bed_cols)
    selectInput(
      "extra_download_cols",
      "Optional extra columns to include",
      choices = extra_choices,
      multiple = TRUE,
      selected = intersect("gene", extra_choices)
    )
  })

  markers_download_table <- reactive({
    df <- filtered_markers()
    sel <- input$marker_table_rows_selected
    if (is.null(sel) || length(sel) == 0) {
      return(tibble(
        `#chrom` = character(0),
        start = integer(0),
        end = integer(0),
        gene_id = character(0),
        aa_position = integer(0),
        strand = character(0)
      ))
    }
    selected_df <- df[sel, , drop = FALSE]
    extra <- input$extra_download_cols %||% character(0)
    extra <- intersect(extra, names(selected_df))
    selected_df %>%
      mutate(
        `#chrom` = .data$chrom,
        end = .data$stop
      ) %>%
      select(`#chrom`, start, end, gene_id, aa_position, strand, all_of(extra))
  })

  output$markers_download_preview <- renderDT({
    datatable(
      markers_download_table(),
      width = "100%",
      options = list(
        pageLength = 10,
        autoWidth = FALSE,
        scrollX = TRUE,
        scrollCollapse = TRUE
      )
    )
  })
  output$download_markers <- downloadHandler(
    filename = function() {
      "selected_markers.bed"
    },
    content = function(file) {
      write_tsv(markers_download_table(), file)
    },
    contentType = "text/plain"
  )

  # Markers selected in the Markers tab — used as the pool for grouping
  selected_markers_for_grouping <- reactive({
    df <- filtered_markers()
    sel <- input$marker_table_rows_selected
    if (is.null(sel) || length(sel) == 0) {
      return(df[0, c("gene_id", "gene", "aa_position"), drop = FALSE])
    }
    df[sel, c("gene_id", "gene", "aa_position"), drop = FALSE]
  })

  output$grouping_placeholder_ui <- renderUI({
    sm <- selected_markers_for_grouping()
    if (nrow(sm) == 0) {
      div(
        class = "alert alert-info",
        role = "alert",
        strong("No markers selected."),
        " Go to the ", strong("Markers"), " tab, select markers, then return here to organize them into groups."
      )
    }
  })

  output$grouping_marker_table <- renderDT({
    sm <- selected_markers_for_grouping()
    datatable(
      sm,
      selection = "multiple",
      width = "100%",
      options = list(
        pageLength = 25,
        autoWidth = FALSE,
        scrollX = TRUE,
        scrollCollapse = TRUE
      )
    )
  })

  # Marker groups: group_id, gene_id, aa_position
  groups <- reactiveVal(tibble(
    group_id = character(),
    gene_id = character(),
    aa_position = numeric()
  ))

  output$groups_download_cols_ui <- renderUI({
    df <- annotated_markers()
    base_cols <- c("gene_id", "aa_position")
    extra_choices <- setdiff(names(df), base_cols)
    selectInput(
      "extra_groups_download_cols",
      "Optional extra columns to include in groups file",
      choices = extra_choices,
      multiple = TRUE
    )
  })

  groups_download_table <- reactive({
    g <- groups()
    if (nrow(g) == 0) {
      return(g)
    }
    df <- annotated_markers()
    extra <- input$extra_groups_download_cols %||% character(0)
    extra <- intersect(extra, names(df))
    df_keyed <- df %>% mutate(gene_id = as.character(.data$gene_id))
    g %>%
      left_join(df_keyed, by = c("gene_id", "aa_position")) %>%
      select(group_id, gene_id, aa_position, all_of(extra))
  })

  observeEvent(input$autogenerate_group_name, {
    if (!isTRUE(input$autogenerate_group_name)) {
      return()
    }
    sm <- selected_markers_for_grouping()
    sel <- input$grouping_marker_table_rows_selected
    if (nrow(sm) == 0 || is.null(sel) || length(sel) == 0) {
      updateCheckboxInput(session, "autogenerate_group_name", value = FALSE)
      showModal(modalDialog(
        title = "Cannot generate a name",
        "Select at least one marker in the table above.",
        easyClose = TRUE
      ))
    }
  }, ignoreInit = TRUE)

  observe({
    if (!isTRUE(input$autogenerate_group_name)) {
      return()
    }
    sm <- selected_markers_for_grouping()
    sel <- input$grouping_marker_table_rows_selected
    if (is.null(sel) || length(sel) == 0 || nrow(sm) == 0) {
      return()
    }
    updateTextInput(
      session,
      "group_name",
      value = stave_group_name(sm[sel, , drop = FALSE])
    )
  })

  observeEvent(input$add_to_group, {
    sm <- selected_markers_for_grouping()
    sel <- input$grouping_marker_table_rows_selected
    group_name <- trimws(input$group_name)
    if (nrow(sm) == 0 || is.null(sel) || length(sel) == 0 || group_name == "") {
      showModal(modalDialog(
        title = "Cannot add to group",
        "Please select at least one marker from the table above and provide a group name.",
        easyClose = TRUE
      ))
      return()
    }
    new_rows <- sm[sel, c("gene_id", "aa_position"), drop = FALSE]
    new_rows <- new_rows %>%
      mutate(
        group_id = group_name,
        gene_id = as.character(.data$gene_id)
      ) %>%
      select(group_id, gene_id, aa_position)
    groups(bind_rows(groups(), new_rows))
  })

  observeEvent(input$clear_groups, {
    groups(tibble(
      group_id = character(),
      gene_id = character(),
      aa_position = numeric()
    ))
  })

  output$group_table <- renderDT({
    datatable(
      groups_download_table(),
      width = "100%",
      options = list(
        pageLength = 25,
        autoWidth = FALSE,
        scrollX = TRUE,
        scrollCollapse = TRUE
      )
    )
  })

  output$download_groups <- downloadHandler(
    filename = function() {
      "marker_groups.tsv"
    },
    content = function(file) {
      write_tsv(groups_download_table(), file)
    }
  )
}

shinyApp(ui, server)
