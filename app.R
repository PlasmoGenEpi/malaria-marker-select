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
      h1("Malaria Resistance Marker Selector"),
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
        helpText("Narrow which markers appear in the table."),
        h4("WHO"),
        helpText(
          "Limit markers using WHO compendium annotations.",
          "Details are on the Instructions tab."
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
          "Upload a BED (chrom, start, end) to keep markers that overlap your panel."
        ),
        fileInput(
          "panel_bed_file",
          "Upload panel BED",
          accept = c(".bed", ".tsv", ".bed.tsv", "text/tab-separated-values", "text/plain")
        ),
        actionButton("clear_panel_bed", "Clear panel BED"),
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
        tags$summary("Custom data"),
        helpText("Optional — replace the built-in marker table with your own TSV."),
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
            "Select individual markers to download or group. Use the sidebar to filter first if needed."
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
            "Exports ", code("#chrom"), ", ", code("start"), ", ", code("end"), ", ",
            code("gene_id"), ", ", code("strand"), ", plus any optional columns you pick."
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
            "Next: open ", strong("Marker groups"), " to organize this selection."
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
            "Select rows above, name the group, then add. Markers can be in more than one group."
          ),
          textInput("group_name", "Group name", value = "", placeholder = "e.g. pfdhfr_pfdhps"),
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
          h3("Overview"),
          p(
            "Browse Plasmodium falciparum drug‑resistance markers, filter by the WHO",
            "Compendium or your panel BED, select markers, download a BED‑like TSV,",
            "and optionally assign them to named groups."
          ),
          p(
            strong("Workflow:"),
            " sidebar filters → select rows on Markers → download (optional) → Marker groups."
          ),
          h3("Using the app"),
          tags$ol(
            tags$li("Use the sidebar to filter by WHO drug/classification and/or panel BED overlap."),
            tags$li("On ", strong("Markers"), ", search or sort the table, then select rows (or Select all)."),
            tags$li("Download selected markers as a BED‑like TSV (", code("#chrom"), ", ", code("start"), ", ", code("end"), ", ", code("gene_id"), ", ", code("strand"), ", plus optional columns)."),
            tags$li("On ", strong("Marker groups"), ", assign selected markers to group names and download a groups TSV (", code("group_id"), ", ", code("gene_id"), ", ", code("aa_position"), ").")
          ),
          h3("Filters"),
          tags$ul(
            tags$li(
              strong("WHO:"),
              " keep markers in the compendium, and/or those linked to selected drugs or",
              " classifications (Validated / Candidate / Potential). A site can appear under",
              " multiple drugs or classes; filters keep it if ", em("any"), " matching WHO entry exists."
            ),
            tags$li(
              strong("Panel BED:"),
              " upload ", code("chrom"), "/", code("start"), "/", code("end"),
              " (header optional). Keeps markers overlapping any interval. Chromosome names must",
              " match the marker table (e.g. ", code("Pf3D7_13_v3"), ")."
            )
          ),
          h3("WHO field guide"),
          p(
            "From the WHO Compendium user guide (Dec 2025, v1.0). Classifications reflect",
            " laboratory, clinical, and/or genetic epidemiological evidence."
          ),
          tags$ul(
            tags$li(strong("Gene(s):"), " gene(s) carrying the alteration."),
            tags$li(strong("Alteration(s):"), " mutation, haplotype, or amplification (", code("who_alterations"), ")."),
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
          ),
          h3("Custom data"),
          p(
            "Optional: upload a marker TSV with ",
            code("gene_id"), ", ", code("gene"), ", ", code("aa_position"), ", ",
            code("chrom"), ", ", code("start"), ", ", code("stop"), ", ", code("strand"),
            ". WHO columns join automatically when gene + aa_position match."
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
    base_bed_cols <- c("chrom", "start", "stop", "strand", "gene_id")
    extra_choices <- setdiff(names(df), base_bed_cols)
    selectInput(
      "extra_download_cols",
      "Optional extra columns to include",
      choices = extra_choices,
      multiple = TRUE,
      selected = intersect(who_display_cols, extra_choices)
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
      select(`#chrom`, start, end, gene_id, strand, all_of(extra))
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
      "selected_markers.bed.tsv"
    },
    content = function(file) {
      write_tsv(markers_download_table(), file)
    }
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
        " Go to the ", strong("Markers"), " tab, select markers (e.g. for download), then return here to organize them into groups."
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
