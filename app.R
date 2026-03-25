library(shiny)
library(DT)
library(readr)
library(dplyr)

# Helper to load the default marker table bundled with the app
load_default_markers <- function() {
  path <- file.path("data", "principal_resistance_marker_info_table.tsv")
  read_tsv(path, show_col_types = FALSE, progress = FALSE)
}

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      .app-header {
        display: flex;
        align-items: center;
        gap: 12px;
      }
      .app-header img {
        height: 60px;
      }
    "))
  ),
  div(
    class = "app-header",
    tags$img(src = "logo.png", alt = "PlasmoGenEpi logo"),
    h1("Malaria Resistance Marker Selector")
  ),
  sidebarLayout(
    sidebarPanel(
      h4("Data source"),
      helpText(
        "This app starts with a built‑in table of principal resistance markers.",
        "You can optionally upload your own table in TSV (tab‑separated) format."
      ),
      fileInput(
        "marker_file",
        "Upload custom marker table (TSV)",
        accept = c(".tsv", "text/tab-separated-values")
      ),
      actionButton("reset_data", "Use built‑in marker table"),
      tags$hr(),
      h4("Table Display"),
      helpText("By default you see gene_id, gene (gene name), and amino acid position."),
      uiOutput("extra_cols_ui"),
      width = 3
    ),
    mainPanel(
      tabsetPanel(
        tabPanel(
          "Markers",
          h3("Browse and select markers"),
          helpText(
            "1. Use the search boxes to find markers of interest (for example by gene).",
            tags$br(),
            "2. Click on rows to select markers; use Ctrl / Cmd or Shift for multiple selections."
          ),
          DTOutput("marker_table"),
          tags$hr(),
          h3("Download selected markers as BED‑like table"),
          helpText(
            "The download always includes #chrom, start, end, gene_id, and strand.",
            "You can optionally add any other columns from the marker table."
          ),
          uiOutput("download_cols_ui"),
          checkboxInput("show_markers_preview", "Show preview before downloading", value = FALSE),
          conditionalPanel(
            condition = "input.show_markers_preview == true",
            h4("Preview of downloaded marker table"),
            helpText("This preview shows exactly what will be written to the downloaded file, based on your current selection and chosen columns."),
            DTOutput("markers_download_preview")
          ),
          downloadButton("download_markers", "Download selected markers"),
          helpText(
            strong("Next step:"), "Go to the ", strong("Marker groups"), " tab to organize your selected markers into groups."
          )
        ),
        tabPanel(
          "Marker groups",
          h3("Organize selected markers into groups"),
          helpText(
            "Markers shown below are the ones you selected in the ", strong("Markers"), " tab.",
            "If the table is empty, go back to the Markers tab, select markers, then return here."
          ),
          uiOutput("grouping_placeholder_ui"),
          DTOutput("grouping_marker_table"),
          tags$hr(),
          h3("Create groups"),
          helpText(
            "Select one or more markers from the table above, enter a group name, then click Add.",
            "You can create as many groups as you like. A marker can belong to multiple groups.",
            "You do not need to assign every selected marker to a group."
          ),
          textInput("group_name", "Group name", value = "", placeholder = "e.g. pfdhfr_pfdhps"),
          actionButton("add_to_group", "Add selected markers to group"),
          actionButton("clear_groups", "Clear all groups"),
          tags$br(), tags$br(),
          h4("Preview of marker-groups file"),
          helpText(
            "This preview shows exactly what will be written to the downloaded marker-groups TSV.",
            "You can optionally include extra columns from the marker table."
          ),
          uiOutput("groups_download_cols_ui"),
          DTOutput("group_table"),
          downloadButton("download_groups", "Download marker groups")
        ),
        tabPanel(
          "Instructions",
          h3("What this app does"),
          p(
            "This app helps you explore a curated list of Plasmodium falciparum drug‑resistance markers,",
            "select markers of interest, download them as a BED‑style file, and optionally organize them into groups."
          ),
          p(
            strong("Workflow:"), " Markers tab → select markers → download (optional) → Marker groups tab → create groups from your selection."
          ),
          h3("Input data format"),
          p(
            "The built‑in marker table includes at least the following columns:",
            code("gene_id"), ", ", code("gene"), ", ", code("aa_position"), ", ",
            code("chrom"), ", ", code("start"), ", ", code("stop"), ", ", code("strand"), "."
          ),
          p(
            "If you upload your own table it should also be a tab‑separated (TSV) file with the same",
            "minimum columns. Any extra columns will be available to display and to include in downloads."
          ),
          h3("Viewing and filtering markers"),
          tags$ul(
            tags$li("Go to the ", strong("Markers"), " tab."),
            tags$li("Type into the search boxes below each column header to filter by gene, position or other fields."),
            tags$li("Click on column headers to sort.")
          ),
          h3("Selecting markers"),
          tags$ul(
            tags$li("Click on a row to select a marker."),
            tags$li("Hold Ctrl (Windows/Linux) or Command (macOS) and click to add or remove individual markers from the selection."),
            tags$li("Use Shift‑click to select a block of markers.")
          ),
          h3("Downloading markers as BED‑like file"),
          tags$ul(
            tags$li(
              "Under ", strong("Download selected markers as BED‑like table"), ", choose any extra columns you would like to include."
            ),
            tags$li(
              "Click ", strong("Download selected markers"), ". The file will contain at least:",
              code("#chrom"), ", ", code("start"), ", ", code("end"), ", ", code("gene_id"), ", ", code("strand"), "."
            )
          ),
          h3("Creating marker groups"),
          tags$ol(
            tags$li("In the ", strong("Markers"), " tab, select the markers you are interested in (and optionally download them as a BED file)."),
            tags$li("Go to the ", strong("Marker groups"), " tab."),
            tags$li("You will see the markers you selected. Choose which of these you want to add to a group."),
            tags$li("Enter a group name (for example ", code("pfdhfr_pfdhps"), ") and click ", strong("Add selected markers to group"), "."),
            tags$li("Repeat to create as many groups as you like. Each marker can appear in multiple groups. You do not need to assign every selected marker to a group.")
          ),
          p(
            "The groups table shows the mapping between groups and individual markers.",
            "You can download this table for downstream analyses."
          ),
          h3("Groups output format"),
          p(
            "The marker‑group file you download is a TSV with these core columns: ",
            code("group_id"), ", ", code("gene_id"), ", ", code("aa_position"), ". ",
            "Optionally, you can include additional columns from your marker table (for example ", code("chrom"), " or ", code("strand"), "). ",
            "This is compatible with the style of the example ", code("example_loci_groups.tsv"), " file."
          ),
          h3("Running and sharing the app"),
          tags$ol(
            tags$li("Install R (version 4 or later) and RStudio if you prefer a graphical interface."),
            tags$li("Install the R packages: ", code('install.packages(c("shiny", "DT", "readr", "dplyr"))'), "."),
            tags$li("Open the project folder in R or RStudio and run: ", code('shiny::runApp()'), "."),
            tags$li(
              "To share online, you can deploy this folder to ", a("shinyapps.io", href = "https://www.shinyapps.io", target = "_blank"),
              " using the ", code("rsconnect"), " package."
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

  # Default and extra columns for display
  base_display_cols <- c("gene_id", "gene", "aa_position")

  output$extra_cols_ui <- renderUI({
    df <- markers()
    all_cols <- names(df)
    extra_choices <- setdiff(all_cols, base_display_cols)
    checkboxGroupInput(
      "extra_cols",
      "Select additional columns to display",
      choices = extra_choices,
      selected = character(0)
    )
  })

  output$marker_table <- renderDT({
    df <- markers()
    display_cols <- unique(c(base_display_cols, input$extra_cols %||% character(0)))
    display_cols <- intersect(display_cols, names(df))
    datatable(
      df[, display_cols, drop = FALSE],
      selection = "multiple",
      filter = "top",
      options = list(pageLength = 25, autoWidth = TRUE)
    )
  })

  # Extra columns for BED-like download
  output$download_cols_ui <- renderUI({
    df <- markers()
    base_bed_cols <- c("chrom", "start", "stop", "strand", "gene_id")
    extra_choices <- setdiff(names(df), base_bed_cols)
    selectInput(
      "extra_download_cols",
      "Optional extra columns to include",
      choices = extra_choices,
      multiple = TRUE
    )
  })

  markers_download_table <- reactive({
    df <- markers()
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
      options = list(pageLength = 10, autoWidth = TRUE)
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
    df <- markers()
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
      options = list(pageLength = 25, autoWidth = TRUE)
    )
  })

  # Marker groups: group_id, gene_id, aa_position
  groups <- reactiveVal(tibble(
    group_id = character(),
    gene_id = character(),
    aa_position = numeric()
  ))

  output$groups_download_cols_ui <- renderUI({
    df <- markers()
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
    df <- markers()
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
      options = list(pageLength = 25, autoWidth = TRUE)
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

