# ==============================================================================
# R SHINY APPLICATION FOR SINGLE-CELL LR INTERACTIONS & PATHWAY VISUALIZATIONS
# ==============================================================================
# Styled with a Premium Modern Clean-Tech Theme (bslib + custom CSS)
# Integrates SingleCellSignalR (NoNet), Reactome Hierarchies, and Day Comparisons.
# ==============================================================================

# 1. Package Checks & Loading

library("shiny")
library("bslib")
library("DT")
library("SingleCellSignalR")
library("BulkSignalR")
library("Seurat")
library("ggplot2")
library("plotly")
library("ggraph")
library("tidygraph")
library("tidyverse")
library("glue")
library("stringr")
library("dplyr")
library("purrr")
library("viridis")
library("RColorBrewer")
library("htmlwidgets")
library("treemap")
library("d3treeR")
library("data.tree")
library("circlepackeR")
library("circlize")

# Increase upload size to 2GB for RDS objects
options(shiny.maxRequestSize = 2000 * 1024^2)

# ==============================================================================
# 2. HELPER FUNCTIONS
# ==============================================================================

# Extract matrix from Seurat object robustly across v3, v4, and v5
get_seurat_matrix <- function(seurat_obj, layer_name = "log_counts", assay_name = "RNA") {
  mat <- NULL

  # Try Seurat v5 layer format first
  try(
    {
      mat <- Seurat::GetAssayData(seurat_obj, assay = assay_name, layer = layer_name)
    },
    silent = TRUE
  )

  # Try Seurat v4 slot format if v5 failed
  if (is.null(mat) || length(mat) == 0) {
    try(
      {
        slot_name <- if (layer_name == "log_counts") "data" else "counts"
        mat <- Seurat::GetAssayData(seurat_obj, assay = assay_name, slot = slot_name)
      },
      silent = TRUE
    )
  }

  # Try Default Assay and "data" slot as fallback
  if (is.null(mat) || length(mat) == 0) {
    try(
      {
        default_assay <- Seurat::DefaultAssay(seurat_obj)
        mat <- Seurat::GetAssayData(seurat_obj, assay = default_assay, slot = "data")
      },
      silent = TRUE
    )
  }

  # Try "counts" slot if "data" is not populated
  if (is.null(mat) || length(mat) == 0) {
    try(
      {
        mat <- Seurat::GetAssayData(seurat_obj, slot = "counts")
      },
      silent = TRUE
    )
  }

  if (is.null(mat) || length(mat) == 0) {
    stop("Could not extract expression matrix from the Seurat object.")
  }

  mat <- as.matrix(mat)
  rownames(mat) <- toupper(rownames(mat))
  return(mat)
}

# Highly optimized, vectorized Reactome pathway mapping
match_common_pathways_vectorized <- function(df_interactions, reactome_human) {
  if (nrow(df_interactions) == 0) {
    return(df_interactions |> dplyr::mutate(
      pathway_name = character(), root_parent = character(),
      subparent = character(), subsubparent = character()
    ))
  }

  # Columns to retain
  reactome_L <- reactome_human |>
    dplyr::select(gene_name, pathway_name, root_parent, subparent, subsubparent) |>
    dplyr::rename(L = gene_name, pathway_name_L = pathway_name)

  reactome_R <- reactome_human |>
    dplyr::select(gene_name, pathway_name, root_parent, subparent, subsubparent) |>
    dplyr::rename(R = gene_name, pathway_name_R = pathway_name)

  # Merge interactions with ligand pathways
  m1 <- df_interactions |>
    dplyr::select(L, R) |>
    dplyr::distinct() |>
    dplyr::inner_join(reactome_L, by = "L", relationship = "many-to-many")

  # Merge with receptor pathways on shared pathway name and parent hierarchies
  m2 <- m1 |>
    dplyr::inner_join(reactome_R, by = c("R", "root_parent", "subparent", "subsubparent"), relationship = "many-to-many") |>
    dplyr::filter(pathway_name_L == pathway_name_R) |>
    dplyr::rename(pathway_name = pathway_name_L) |>
    dplyr::select(L, R, pathway_name, root_parent, subparent, subsubparent)

  if (nrow(m2) > 0) {
    # Collapse multiple matching common pathways with a semicolon
    m2_collapsed <- m2 |>
      dplyr::group_by(L, R) |>
      dplyr::summarise(
        pathway_name = paste(pathway_name, collapse = "; "),
        root_parent = paste(root_parent, collapse = "; "),
        subparent = paste(subparent, collapse = "; "),
        subsubparent = paste(subsubparent, collapse = "; "),
        .groups = "drop"
      )

    # Left join collapsed paths back to original interactions
    res <- df_interactions |>
      dplyr::left_join(m2_collapsed, by = c("L", "R"))
  } else {
    # No pathways matched
    res <- df_interactions |>
      dplyr::mutate(
        pathway_name = NA_character_,
        root_parent = NA_character_,
        subparent = NA_character_,
        subsubparent = NA_character_
      )
  }

  return(res)
}

# Safe wrapper around getParacrines/getAutocrines
get_pair_interactions_safe <- function(scsrnn, sender, receiver, reactome_human) {
  inter <- NULL
  tryCatch(
    {
      if (sender == receiver) {
        inter <- getAutocrines(scsrnn, sender)
      } else {
        inter <- getParacrines(scsrnn, sender, receiver)
      }
    },
    error = function(e) {
      # If not in inferences, return NULL
      inter <- NULL
    }
  )

  if (!is.null(inter) && nrow(inter) > 0) {
    inter$Sender <- sender
    inter$Receiver <- receiver

    # Map pathways using our vectorized logic
    inter <- match_common_pathways_vectorized(inter, reactome_human)
    return(inter)
  }
  return(data.frame())
}

# Helper to compile multiple cell-cell pairs dynamically
get_multiple_interactions_safe <- function(scsrnn, senders, receivers, reactome_human) {
  results_list <- list()
  for (s in senders) {
    for (r in receivers) {
      df_pair <- get_pair_interactions_safe(scsrnn, s, r, reactome_human)
      if (nrow(df_pair) > 0) {
        results_list[[length(results_list) + 1]] <- df_pair
      }
    }
  }
  if (length(results_list) > 0) {
    return(bind_rows(results_list))
  }
  return(data.frame())
}

# Helper to compile all interactions for a given day
get_all_interactions_for_day <- function(scsrnn, clusters, reactome_human) {
  results_list <- list()
  for (c1 in clusters) {
    for (c2 in clusters) {
      tryCatch(
        {
          if (c1 == c2) {
            inter <- getAutocrines(scsrnn, c1)
          } else {
            inter <- getParacrines(scsrnn, c1, c2)
          }
          if (!is.null(inter) && nrow(inter) > 0) {
            inter$Sender <- c1
            inter$Receiver <- c2
            results_list[[length(results_list) + 1]] <- inter
          }
        },
        error = function(e) {
          # ignore missing DE pairs
        }
      )
    }
  }
  if (length(results_list) > 0) {
    combined <- bind_rows(results_list)
    combined <- match_common_pathways_vectorized(combined, reactome_human)
    return(combined)
  }
  return(data.frame())
}

# Helper to retrieve pathway annotations for a specific condition/day and scope
get_pathway_frequency_data <- function(scsrnn, scope, sender, receiver, reactome_human,
                                       min_lr_score = 0.7, max_pval = 0.05, min_logFC = 0.05) {
  if (scope == "all") {
    clusters <- unique(as.character(slot(scsrnn, "populations")))
    df_pair <- get_all_interactions_for_day(scsrnn, clusters, reactome_human)
  } else {
    if (is.null(sender) || is.null(receiver)) {
      return(data.frame())
    }
    df_pair <- get_pair_interactions_safe(scsrnn, sender, receiver, reactome_human)
  }

  if (is.null(df_pair) || nrow(df_pair) == 0) {
    return(data.frame())
  }

  # Filter by thresholds
  df_pair <- df_pair |>
    dplyr::filter(LR.score >= min_lr_score) |>
    dplyr::filter(pval <= max_pval) |>
    dplyr::filter(L.logFC >= min_logFC & R.logFC >= min_logFC)

  if (nrow(df_pair) == 0) {
    return(data.frame())
  }

  concat_csv_subparent_root <- df_pair |>
    dplyr::mutate(pathway_name = dplyr::if_else(pathway_name == "" | is.na(pathway_name), NA_character_, pathway_name)) |>
    tidyr::drop_na(pathway_name) |>
    dplyr::select(root_parent, subparent) |>
    dplyr::mutate(
      root_parent = stringr::str_split(root_parent, "; "),
      subparent = stringr::str_split(subparent, "; ")
    ) |>
    tidyr::unnest(c(root_parent, subparent))

  return(concat_csv_subparent_root)
}

# ==============================================================================
# 3. SHINY USER INTERFACE (UI)
# ==============================================================================

# Custom Premium Styling CSS
custom_css <- "
  /* High-end CSS for modern premium dashboard */
  body {
    font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
    background-color: #f8fafc;
    color: #0f172a;
  }
  .navbar {
    background: linear-gradient(135deg, #0f172a 0%, #1e1b4b 50%, #312e81 100%) !important;
    box-shadow: 0 4px 20px -2px rgba(0, 0, 0, 0.15);
    border-bottom: 1px solid rgba(255, 255, 255, 0.08);
  }
  .navbar-brand {
    font-weight: 800 !important;
    letter-spacing: -0.7px;
    background: linear-gradient(to right, #a5b4fc, #38bdf8);
    -webkit-background-clip: text;
    -webkit-text-fill-color: transparent;
  }
  .navbar .nav-link {
    font-weight: 600;
    color: #cbd5e1 !important;
    transition: all 0.25s ease;
    border-radius: 6px;
    margin: 0 3px;
  }
  .navbar .nav-link:hover, .navbar .nav-link.active {
    color: #ffffff !important;
    background-color: rgba(255, 255, 255, 0.08);
    text-shadow: 0 0 10px rgba(165, 180, 252, 0.4);
  }
  .card {
    border-radius: 14px !important;
    border: 1px solid #e2e8f0 !important;
    box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.04), 0 2px 4px -1px rgba(0, 0, 0, 0.02) !important;
    transition: transform 0.2s ease, box-shadow 0.2s ease !important;
    background: rgba(255, 255, 255, 0.95);
    backdrop-filter: blur(10px);
  }
  .card:hover {
    box-shadow: 0 12px 20px -3px rgba(0, 0, 0, 0.07), 0 4px 8px -2px rgba(0, 0, 0, 0.03) !important;
  }
  .card-header {
    background-color: #f8fafc !important;
    border-bottom: 1px solid #e2e8f0 !important;
    font-weight: 700 !important;
    color: #1e293b !important;
    font-size: 1.05rem;
    padding: 12px 20px !important;
    border-top-left-radius: 14px !important;
    border-top-right-radius: 14px !important;
  }
  .btn-primary {
    background: linear-gradient(135deg, #4f46e5 0%, #6366f1 100%) !important;
    border: none !important;
    font-weight: 600 !important;
    letter-spacing: -0.2px;
    padding: 10px 20px !important;
    border-radius: 8px !important;
    box-shadow: 0 4px 12px -1px rgba(79, 70, 229, 0.35) !important;
    transition: all 0.2s ease !important;
  }
  .btn-primary:hover {
    background: linear-gradient(135deg, #4338ca 0%, #4f46e5 100%) !important;
    transform: translateY(-1.5px) !important;
    box-shadow: 0 6px 16px -1px rgba(79, 70, 229, 0.45) !important;
  }
  .btn-success {
    background: linear-gradient(135deg, #059669 0%, #10b981 100%) !important;
    border: none !important;
    font-weight: 600 !important;
    border-radius: 8px !important;
    box-shadow: 0 4px 12px -1px rgba(5, 150, 105, 0.35) !important;
  }
  .btn-success:hover {
    background: linear-gradient(135deg, #047857 0%, #059669 100%) !important;
    transform: translateY(-1.5px) !important;
  }
  .shiny-notification {
    background-color: #1e1b4b !important;
    border: 1px solid #4f46e5 !important;
    color: #ffffff !important;
    border-radius: 10px !important;
    opacity: 0.95 !important;
  }
  .sidebar {
    background-color: #ffffff !important;
    border-right: 1px solid #e2e8f0 !important;
    padding: 24px !important;
  }
  .value-box {
    border-radius: 12px;
    border: 1px solid #e2e8f0;
    padding: 16px;
    background: #ffffff;
    box-shadow: 0 2px 4px rgba(0, 0, 0, 0.02);
  }
  .badge-count {
    background-color: #e0e7ff;
    color: #4338ca;
    padding: 4px 8px;
    border-radius: 9999px;
    font-size: 0.8rem;
    font-weight: 600;
  }
  .d3tree2 rect {
    stroke: black !important;
    stroke-width: 1px !important;
  }
  .node--leaf {
    fill: inherit !important;
    fill-opacity: 0.2 !important;
    stroke: rgba(0, 0, 0, 1) !important;
    stroke-width: 0.5px !important;
  }
"


ui <- do.call(page_navbar, c(
  list(
    title = "SingleCellSignalR Portal",
    theme = bs_theme(
      version = 5,
      bootswatch = "flatly",
      primary = "#4f46e5",
      secondary = "#06b6d4"
    ),
    header = tags$head(
      tags$style(HTML(custom_css)),
      tags$link(rel = "stylesheet", href = "https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700;800&display=swap"),
      tags$script(HTML("
        function downloadWidgetAsPNG(widgetId, filename) {
          var el = document.getElementById(widgetId);
          if (!el) return;
          var svg = el.querySelector('svg');
          if (!svg) {
            alert('No SVG element found in this widget.');
            return;
          }
          var bbox = svg.getBoundingClientRect();
          var width = svg.getAttribute('width') || bbox.width || 800;
          var height = svg.getAttribute('height') || bbox.height || 800;
          
          var svgClone = svg.cloneNode(true);
          svgClone.setAttribute('width', width);
          svgClone.setAttribute('height', height);
          svgClone.style.backgroundColor = '#ffffff';
          
          var svgString = new XMLSerializer().serializeToString(svgClone);
          var svgBlob = new Blob([svgString], {type: 'image/svg+xml;charset=utf-8'});
          var URL = window.URL || window.webkitURL || window;
          var blobURL = URL.createObjectURL(svgBlob);
          
          var image = new Image();
          image.onload = function() {
            var canvas = document.createElement('canvas');
            canvas.width = width * 2;
            canvas.height = height * 2;
            var context = canvas.getContext('2d');
            context.scale(2, 2);
            context.fillStyle = '#ffffff';
            context.fillRect(0, 0, width, height);
            context.drawImage(image, 0, 0, width, height);
            
            var pngURL = canvas.toDataURL('image/png');
            var downloadLink = document.createElement('a');
            downloadLink.href = pngURL;
            downloadLink.download = filename;
            document.body.appendChild(downloadLink);
            downloadLink.click();
            document.body.removeChild(downloadLink);
            URL.revokeObjectURL(blobURL);
          };
          image.src = blobURL;
        }
      "))
    )
  ),
  list(
    # ========================================================================
    # TAB 1: DATA LOADING & RUN CONFIGURATION
    # ========================================================================
    nav_panel(
      "Data Loading & Config",
      layout_sidebar(
        sidebar = sidebar(
          title = "Configuration Panel",
          width = 360,

          # DataSource Select
          card(
            card_header("Data Source Selection"),
            radioButtons("data_source", "Select Source:",
              choices = c(
                "Demo Workspace scsrnn Objects" = "demo",
                "Upload scsrnn RDS Files" = "upload"
              ),
              selected = "demo"
            ),
            conditionalPanel(
              condition = "input.data_source == 'upload'",
              fileInput("scsrnn_files", "Upload scsrnn RDS File(s):", multiple = TRUE, accept = c(".rds"))
            )
          ),

          # Inference Parameters
          card(
            card_header("Filtering Thresholds"),
            sliderInput("min_LR_score", "Min LR Score (Post-filter):", min = 0.7, max = 1, value = 0.7, step = 0.05),
            sliderInput("max_pval", "Max P-Value threshold:", min = 0, max = 0.05, value = 0.05, step = 0.01),
            numericInput("min_logFC", "Min Log Fold Change:", value = 0.05, min = 0, step = 0.01)
          ),

          # Run Action
          actionButton("run_pipeline", "Load and Initialize Datasets", class = "btn-primary w-100 mt-2")
        ),

        # Main Content Data loading
        layout_column_wrap(
          width = 1,
          card(
            card_header("Loaded scsrnn Objects Information"),
            uiOutput("dataset_stats"),
            hr(),
            h5("Pipeline Execution Log"),
            verbatimTextOutput("pipeline_log")
          )
        )
      )
    ),

    # ========================================================================
    # TAB 2: GLOBAL LR INTERACTION EXPLORER
    # ========================================================================
    nav_panel(
      "Global LR Explorer",
      layout_sidebar(
        sidebar = sidebar(
          title = "Interaction Filters",
          width = 360,
          selectInput("explorer_day", "Select Dataset / Day:", choices = NULL),
          selectInput("explorer_sender", "Sender Population:", choices = NULL, multiple = TRUE),
          div(
            style = "margin-top: -12px; margin-bottom: 15px;",
            actionButton("sender_all", "Select All", class = "btn btn-sm btn-outline-secondary", style = "padding: 1px 6px; font-size: 0.7rem; font-weight: 500;"),
            actionButton("sender_none", "Clear All", class = "btn btn-sm btn-outline-secondary", style = "padding: 1px 6px; font-size: 0.7rem; font-weight: 500; margin-left: 5px;")
          ),
          selectInput("explorer_receiver", "Receiver Population:", choices = NULL, multiple = TRUE),
          div(
            style = "margin-top: -12px; margin-bottom: 15px;",
            actionButton("receiver_all", "Select All", class = "btn btn-sm btn-outline-secondary", style = "padding: 1px 6px; font-size: 0.7rem; font-weight: 500;"),
            actionButton("receiver_none", "Clear All", class = "btn btn-sm btn-outline-secondary", style = "padding: 1px 6px; font-size: 0.7rem; font-weight: 500; margin-left: 5px;")
          ),
          card(
            card_header("Filtering Thresholds"),
            sliderInput("explorer_lr_threshold", "LR Score Cutoff:", min = 0.7, max = 1, value = 0.7, step = 0.05),
            sliderInput("explorer_max_pval", "Max P-Value threshold:", min = 0, max = 0.05, value = 0.05, step = 0.01),
            numericInput("explorer_min_logFC", "Min Log Fold Change:", value = 0.05, min = 0, step = 0.01),
            textInput("explorer_gene", "Filter by Gene (HUGO symbol):", placeholder = "e.g. CXCL, IL")
          ),
        ),
        navset_card_tab(
          id = "explorer_tabs",
          nav_panel(
            "Global Frequency Heatmap",
            plotOutput("global_heatmap", height = "550px"),
            div(style = "text-align: right; margin-top: 10px; margin-right: 10px; margin-bottom: 5px;",
                downloadButton("download_heatmap_png", "Download Plot (PNG)", class = "btn btn-sm btn-outline-primary"))
          ),
          nav_panel(
            "Network Bubble Plot (Populations)",
            plotOutput("network_bubble", height = "550px"),
            div(style = "text-align: right; margin-top: 10px; margin-right: 10px; margin-bottom: 5px;",
                downloadButton("download_network_bubble_png", "Download Plot (PNG)", class = "btn btn-sm btn-outline-primary"))
          ),
          nav_panel(
            "Detailed L-R Bubble Plot",
            plotOutput("global_bubble", height = "650px"),
            div(style = "text-align: right; margin-top: 10px; margin-right: 10px; margin-bottom: 5px;",
                downloadButton("download_global_bubble_png", "Download Plot (PNG)", class = "btn btn-sm btn-outline-primary"))
          ),
          nav_panel(
            "Interactions Data Table",
            DTOutput("interactions_table")
          )
        )
      )
    ),

    # ========================================================================
    # TAB 3: CHORD DIAGRAM
    # ========================================================================
    nav_panel(
      "Chord Diagram",
      layout_sidebar(
        sidebar = sidebar(
          title = "Chord Diagram Filters",
          width = 360,
          selectInput("chord_day", "Select Dataset / Day:", choices = NULL),
          selectInput("chord_sender", "Sender Population:", choices = NULL, multiple = TRUE),
          div(
            style = "margin-top: -12px; margin-bottom: 15px;",
            actionButton("chord_sender_all", "Select All", class = "btn btn-sm btn-outline-secondary", style = "padding: 1px 6px; font-size: 0.7rem; font-weight: 500;"),
            actionButton("chord_sender_none", "Clear All", class = "btn btn-sm btn-outline-secondary", style = "padding: 1px 6px; font-size: 0.7rem; font-weight: 500; margin-left: 5px;")
          ),
          selectInput("chord_receiver", "Receiver Population:", choices = NULL, multiple = TRUE),
          div(
            style = "margin-top: -12px; margin-bottom: 15px;",
            actionButton("chord_receiver_all", "Select All", class = "btn btn-sm btn-outline-secondary", style = "padding: 1px 6px; font-size: 0.7rem; font-weight: 500;"),
            actionButton("chord_receiver_none", "Clear All", class = "btn btn-sm btn-outline-secondary", style = "padding: 1px 6px; font-size: 0.7rem; font-weight: 500; margin-left: 5px;")
          ),
          card(
            card_header("Filtering Thresholds"),
            sliderInput("chord_lr_threshold", "Min LR Score:", min = 0.7, max = 1, value = 0.7, step = 0.05),
            sliderInput("chord_max_pval", "Max P-Value threshold:", min = 0, max = 0.05, value = 0.05, step = 0.01),
            numericInput("chord_min_logFC", "Min Log Fold Change:", value = 0.05, min = 0, step = 0.01)
          ),
        ),
        card(
          card_header(span(icon("project-diagram"), "Directed Cell-Cell Interaction Chord Diagram")),
          plotOutput("chord_plot", height = "750px"),
          div(style = "text-align: right; margin-top: 10px; margin-right: 15px; margin-bottom: 15px;",
              downloadButton("download_chord_png", "Download Plot (PNG)", class = "btn btn-sm btn-outline-primary"))
        )
      )
    ),

    # ========================================================================
    # TAB 4: PATHWAY DISTRIBUTION
    # ========================================================================
    nav_panel(
      "Pathway Distribution",
      layout_sidebar(
        sidebar = sidebar(
          title = "Pathway Filters",
          width = 360,
          selectInput("dist_day1", "Select Day 1 (Left):", choices = NULL),
          selectInput("dist_day2", "Select Day 2 (Right):", choices = NULL),
          radioButtons("dist_scope", "Analysis Scope:",
            choices = c("Entire Dataset / Day" = "all", "Specific Cell-Cell Pair" = "pair"),
            selected = "all"
          ),
          conditionalPanel(
            condition = "input.dist_scope == 'pair'",
            selectInput("dist_sender", "Sender Cluster:", choices = NULL),
            selectInput("dist_receiver", "Receiver Cluster:", choices = NULL)
          ),
          selectInput("dist_level", "Hierarchy Level:",
            choices = c("Root Parent" = "root_parent", "Subparent" = "subparent"),
            selected = "subparent"
          ),
          card(
            card_header("Filtering Thresholds"),
            sliderInput("dist_lr_threshold", "Min LR Score:", min = 0.7, max = 1, value = 0.7, step = 0.05),
            sliderInput("dist_max_pval", "Max P-Value threshold:", min = 0, max = 0.05, value = 0.05, step = 0.01),
            numericInput("dist_min_logFC", "Min Log Fold Change:", value = 0.05, min = 0, step = 0.01)
          )
        ),
        layout_column_wrap(
          width = 1 / 2,
          card(
            card_header(uiOutput("header_day1")),
            plotlyOutput("pathway_stacked_bar1", height = "650px"),
            div(style = "text-align: right; margin-top: 10px; margin-right: 15px; margin-bottom: 15px;",
                downloadButton("download_stacked1_png", "Download Plot (PNG)", class = "btn btn-sm btn-outline-primary"))
          ),
          card(
            card_header(uiOutput("header_day2")),
            plotlyOutput("pathway_stacked_bar2", height = "650px"),
            div(style = "text-align: right; margin-top: 10px; margin-right: 15px; margin-bottom: 15px;",
                downloadButton("download_stacked2_png", "Download Plot (PNG)", class = "btn btn-sm btn-outline-primary"))
          )
        )
      )
    ),

    # ========================================================================
    # TAB 5: PATHWAY CIRCLEPACK
    # ========================================================================
    nav_panel(
      "Pathway Circlepack",
      layout_sidebar(
        sidebar = sidebar(
          title = "Circlepack Filters",
          width = 360,
          selectInput("circle_day", "Select Dataset / Day:", choices = NULL),
          radioButtons("circle_scope", "Analysis Scope:",
            choices = c("Entire Dataset / Day" = "all", "Specific Cell-Cell Pair" = "pair"),
            selected = "all"
          ),
          conditionalPanel(
            condition = "input.circle_scope == 'pair'",
            selectInput("circle_sender", "Sender Cluster:", choices = NULL),
            selectInput("circle_receiver", "Receiver Cluster:", choices = NULL)
          ),
          card(
            card_header("Filtering Thresholds"),
            sliderInput("circle_lr_threshold", "Min LR Score:", min = 0.7, max = 1, value = 0.7, step = 0.05),
            sliderInput("circle_max_pval", "Max P-Value threshold:", min = 0, max = 0.05, value = 0.05, step = 0.01),
            numericInput("circle_min_logFC", "Min Log Fold Change:", value = 0.05, min = 0, step = 0.01)
          ),
          span("*Note: Circlepack displays Root Parent -> Subparent taxonomy.", style = "font-size: 0.8rem; color: #64748b; font-style: italic;")
        ),
        card(
          card_header(span(icon("dot-circle"), "Circlepack Pathway Hierarchy")),
          circlepackeROutput("pathway_circlepack", height = "750px"),
          div(style = "text-align: right; margin-top: 10px; margin-right: 15px; margin-bottom: 15px;",
              actionButton("download_circle_png", "Download Plot (PNG)", class = "btn btn-sm btn-outline-primary", onclick = "downloadWidgetAsPNG('pathway_circlepack', 'circlepack.png')"))
        )
      )
    ),

    # ========================================================================
    # TAB 6: PATHWAY TREEMAP
    # ========================================================================
    nav_panel(
      "Pathway Treemap",
      layout_sidebar(
        sidebar = sidebar(
          title = "Treemap Filters",
          width = 360,
          selectInput("treemap_day", "Select Dataset / Day:", choices = NULL),
          radioButtons("treemap_scope", "Analysis Scope:",
            choices = c("Entire Dataset / Day" = "all", "Specific Cell-Cell Pair" = "pair"),
            selected = "all"
          ),
          conditionalPanel(
            condition = "input.treemap_scope == 'pair'",
            selectInput("treemap_sender", "Sender Cluster:", choices = NULL),
            selectInput("treemap_receiver", "Receiver Cluster:", choices = NULL)
          ),
          card(
            card_header("Filtering Thresholds"),
            sliderInput("treemap_lr_threshold", "Min LR Score:", min = 0.7, max = 1, value = 0.7, step = 0.05),
            sliderInput("treemap_max_pval", "Max P-Value threshold:", min = 0, max = 0.05, value = 0.05, step = 0.01),
            numericInput("treemap_min_logFC", "Min Log Fold Change:", value = 0.05, min = 0, step = 0.01)
          )
        ),
        card(
          card_header(span(icon("th-large"), "Interactive Treemap (d3treeR)")),
          d3tree2Output("pathway_treemap", height = "750px"),
          div(style = "text-align: right; margin-top: 10px; margin-right: 15px; margin-bottom: 15px;",
              actionButton("download_treemap_png", "Download Plot (PNG)", class = "btn btn-sm btn-outline-primary", onclick = "downloadWidgetAsPNG('pathway_treemap', 'treemap.png')"))
        )
      )
    ),

    # ========================================================================
    # TAB 7: COMPARE UNIQUE INTERACTIONS (2 DAYS)
    # ========================================================================
    nav_panel(
      "Compare Unique Interactions (2 Days)",
      layout_sidebar(
        sidebar = sidebar(
          title = "Comparison Config",
          width = 360,
          selectInput("comp2_day1", "Dataset 1 (Left Side):", choices = NULL),
          selectInput("comp2_day2", "Dataset 2 (Right Side):", choices = NULL),
          selectInput("comp2_sender", "Sender Cluster:", choices = NULL),
          selectInput("comp2_receiver", "Receiver Cluster:", choices = NULL),
          selectInput("comp2_level", "Pathway Hierarchy Level:",
            choices = c("Root Parent" = "root_parent", "Subparent" = "subparent"),
            selected = "root_parent"
          ),
          card(
            card_header("Filtering Thresholds"),
            sliderInput("comp2_lr_threshold", "Min LR Score:", min = 0.7, max = 1, value = 0.7, step = 0.05),
            sliderInput("comp2_max_pval", "Max P-Value threshold:", min = 0, max = 0.05, value = 0.05, step = 0.01),
            numericInput("comp2_min_logFC", "Min Log Fold Change:", value = 0.05, min = 0, step = 0.01)
          )
        ),
        card(
          card_header("Diverging Mirrored Bar Chart (Non-Common Pathways)"),
          p("Compares ligand-receptor interactions unique to each dataset (non-intersecting L-R pairs)."),
          plotOutput("compare2_plot", height = "650px"),
          div(style = "text-align: right; margin-top: 10px; margin-right: 15px; margin-bottom: 15px;",
              downloadButton("download_compare2_png", "Download Plot (PNG)", class = "btn btn-sm btn-outline-primary"))
        )
      )
    ),

    # ========================================================================
    # TAB 8: COMPARE UNIQUE INTERACTIONS (3 DAYS)
    # ========================================================================
    nav_panel(
      "Compare Unique Interactions (3 Days)",
      layout_sidebar(
        sidebar = sidebar(
          title = "Comparison Config",
          width = 360,
          p("Compares pathways across three distinct datasets simultaneously (excluding three-way intersect)."),
          selectInput("comp3_sender", "Sender Cluster:", choices = NULL),
          selectInput("comp3_receiver", "Receiver Cluster:", choices = NULL),
          selectInput("comp3_level", "Pathway Hierarchy Level:",
            choices = c("Root Parent" = "root_parent", "Subparent" = "subparent", "Pathway Name" = "pathway_name"),
            selected = "root_parent"
          ),
          card(
            card_header("Filtering Thresholds"),
            sliderInput("comp3_lr_threshold", "Min LR Score:", min = 0.7, max = 1, value = 0.7, step = 0.05),
            sliderInput("comp3_max_pval", "Max P-Value threshold:", min = 0, max = 0.05, value = 0.05, step = 0.01),
            numericInput("comp3_min_logFC", "Min Log Fold Change:", value = 0.05, min = 0, step = 0.01)
          )
        ),
        card(
          card_header("Pathway Frequency Dumbbell Plot (3-Way Comparison)"),
          plotOutput("compare3_plot", height = "700px"),
          div(style = "text-align: right; margin-top: 10px; margin-right: 15px; margin-bottom: 15px;",
              downloadButton("download_compare3_png", "Download Plot (PNG)", class = "btn btn-sm btn-outline-primary"))
        )
      )
    )
  )
))


# ==============================================================================
# 4. SHINY SERVER LOGIC
# ==============================================================================

server <- function(input, output, session) {
  # Reactive state values
  values <- reactiveValues(
    raw_reactome = NULL, # Raw Reactome reference mapping
    reactome_human = NULL, # End-of-branch filtered Reactome
    scsr_results = list(), # List of inferred SCSRNoNet objects
    interaction_tables = list(), # Collapsed and pathway-mapped interactions
    ref_couleurs = NULL, # Global subparent color lookup dataframe
    max_colors = 1,
    root_colors = NULL, # Global root parent color lookup vector
    pipeline_log = "",
    is_analyzed = FALSE
  )

  # Append progress messages to the UI Log
  log_message <- function(msg) {
    message(msg)
    isolate({
      values$pipeline_log <- paste0(values$pipeline_log, "[", format(Sys.time(), "%H:%M:%S"), "] ", msg, "\n")
    })
  }

  # ----------------------------------------------------------------------------
  # A. LOAD REFERENCE DATABASES ON LAUNCH
  # ----------------------------------------------------------------------------
  observe({
    isolate({
      # 1. Load Reactome Reference
      reactome_path <- "reactome_human_with_gene_names.csv"
      if (file.exists(reactome_path)) {
        log_message("Loading Reactome database reference from workspace...")
        tryCatch(
          {
            values$raw_reactome <- read.csv(reactome_path)
            values$reactome_human <- values$raw_reactome |> dplyr::filter(end_of_branch)
            log_message("Reactome database loaded successfully!")

            # Build global color lookup for all subparents in database!
            log_message("Structuring global pathway color palette reference...")
            ref_couleurs <- values$reactome_human |>
              dplyr::select(root_parent, subparent) |>
              dplyr::distinct() |>
              dplyr::group_by(root_parent) |>
              dplyr::arrange(subparent) |>
              dplyr::mutate(index_couleur = dplyr::row_number()) |>
              dplyr::ungroup()

            values$ref_couleurs <- ref_couleurs
            values$max_colors <- max(ref_couleurs$index_couleur, na.rm = TRUE)
            if (is.infinite(values$max_colors) || is.na(values$max_colors)) values$max_colors <- 1

            # Generate a stable mapping of root_parent to color
            all_roots <- sort(unique(values$reactome_human$root_parent))
            all_roots <- all_roots[all_roots != "" & !is.na(all_roots)]
            root_palette <- colorRampPalette(c(
              RColorBrewer::brewer.pal(12, "Set3"),
              RColorBrewer::brewer.pal(8, "Set2"),
              RColorBrewer::brewer.pal(9, "Set1")
            ))(length(all_roots))
            names(root_palette) <- all_roots
            values$root_colors <- root_palette

            log_message("Pathway color index initialized.")
          },
          error = function(e) {
            log_message(paste("Error loading Reactome CSV:", e$message))
          }
        )
      } else {
        log_message("Warning: reactome_human_with_gene_names.csv not found in app directory!")
        shiny::showNotification("Reactome database file missing. Pathway annotations will be skipped.", type = "warning", duration = NULL)
      }
    })
  })

  # ----------------------------------------------------------------------------
  # B. DISPLAY STATS FOR LOADED scsrnn OBJECTS
  # ----------------------------------------------------------------------------
  output$dataset_stats <- renderUI({
    if (!values$is_analyzed || length(values$scsr_results) == 0) {
      return(div(
        class = "alert alert-warning",
        strong("No Datasets Loaded: "), "Please load the Demo workspace files or upload your scsrnn RDS files."
      ))
    }

    # Create a nice HTML table summarizing each loaded scsrnn object
    rows <- lapply(names(values$scsr_results), function(day) {
      scsrnn <- values$scsr_results[[day]]

      # Extract info safely
      pops <- unique(as.character(slot(scsrnn, "populations")))
      num_pops <- length(pops)
      total_cells <- length(slot(scsrnn, "populations"))

      # Try accessing inferences slot
      paracrine_list <- paracrines(scsrnn)

      num_paracrine_interactions <- if (!is.null(paracrine_list)) {
        sum(sapply(paracrine_list, nrow))
      } else {
        0
      }

      autocrine_list <- autocrines(scsrnn)

      num_autocrine_interactions <- if (!is.null(autocrine_list)) {
        sum(sapply(autocrine_list, nrow))
      } else {
        0
      }

      num_interactions <- num_paracrine_interactions + num_autocrine_interactions

      tags$tr(
        tags$td(strong(day), style = "padding: 10px; border-bottom: 1px solid #e2e8f0;"),
        tags$td(format(total_cells, big.mark = ","), style = "padding: 10px; border-bottom: 1px solid #e2e8f0;"),
        tags$td(format(num_pops, big.mark = ","), style = "padding: 10px; border-bottom: 1px solid #e2e8f0;"),
        tags$td(format(num_interactions, big.mark = ","), style = "padding: 10px; border-bottom: 1px solid #e2e8f0;")
      )
    })

    tags$div(
      class = "table-responsive",
      tags$table(
        class = "table table-hover",
        style = "margin-bottom: 0;",
        tags$thead(
          tags$tr(
            tags$th("Condition / Day", style = "background-color: #f8fafc; padding: 12px; font-weight: 700;"),
            tags$th("Total Cells", style = "background-color: #f8fafc; padding: 12px; font-weight: 700;"),
            tags$th("Number of Populations", style = "background-color: #f8fafc; padding: 12px; font-weight: 700;"),
            tags$th("Total Inferred Interactions", style = "background-color: #f8fafc; padding: 12px; font-weight: 700;")
          )
        ),
        tags$tbody(rows)
      )
    )
  })

  output$pipeline_log <- renderText({
    values$pipeline_log
  })

  # ----------------------------------------------------------------------------
  # C. BACKEND scsrnn DATA LOADING & INITIALIZATION
  # ----------------------------------------------------------------------------
  observeEvent(input$run_pipeline, {
    # 1. Show full-screen blocking loading spinner modal
    showModal(modalDialog(
      title = NULL,
      footer = NULL,
      easyClose = FALSE,
      size = "m",
      div(
        style = "text-align: center; padding: 40px 20px;",
        div(
          class = "spinner-border text-primary",
          role = "status",
          style = "width: 4.5rem; height: 4.5rem; border-width: 0.35em; color: #4f46e5 !important; margin-bottom: 25px;"
        ),
        h3("Running Analysis Pipeline...", style = "font-weight: 800; color: #0f172a; letter-spacing: -0.5px;"),
        p("Calculating cell-cell interactions using SingleCellSignalR NoNet method & mapping Reactome pathway hierarchies.", style = "color: #475569; font-size: 1.05rem; margin-top: 10px;"),
        div(
          class = "progress mt-4",
          style = "height: 8px; border-radius: 4px; background-color: #f1f5f9;",
          div(
            id = "pipeline_progress_bar",
            class = "progress-bar progress-bar-striped progress-bar-animated",
            role = "progressbar",
            style = "width: 100%; background: linear-gradient(to right, #4f46e5, #06b6d4);"
          )
        )
      )
    ))

    # 2. Execute pipeline in a safe tryCatch block
    tryCatch({
      log_message("Starting scsrnn RDS loading process...")
      values$is_analyzed <- FALSE

      scsr_results <- list()

      if (input$data_source == "demo") {
        # Auto-detect workspace scsrnn files
        local_files <- c("scsrnn_step1_local_D3.rds", "scsrnn_step1_local_D7.rds", "scsrnn_step1_local_D14.rds")

        target_files <- c()
        if (all(file.exists(local_files))) {
          target_files <- local_files
        } else {
          # Scan folder for scsrnn rds files
          all_rds <- list.files(pattern = "^scsrnn.*\\.rds$")
          if (length(all_rds) > 0) {
            target_files <- all_rds
          }
        }

        if (length(target_files) == 0) {
          stop("No workspace scsrnn RDS files found. Please upload manually.")
        }

        for (f in target_files) {
          log_message(sprintf("Loading workspace RDS file: %s", f))
          # Extract label
          lbl <- tools::file_path_sans_ext(basename(f))
          lbl <- gsub("^scsrnn_step1_local_", "", lbl)
          lbl <- gsub("^scsrnn_", "", lbl)

          scsr_results[[lbl]] <- readRDS(f)
        }
      } else {
        # Uploaded files
        req(input$scsrnn_files)
        files_df <- input$scsrnn_files

        for (i in seq_len(nrow(files_df))) {
          orig_name <- files_df$name[i]
          path <- files_df$datapath[i]

          log_message(sprintf("Loading uploaded RDS file: %s", orig_name))
          # Extract label
          lbl <- tools::file_path_sans_ext(orig_name)
          lbl <- gsub("^scsrnn_step1_local_", "", lbl)
          lbl <- gsub("^scsrnn_", "", lbl)

          scsr_results[[lbl]] <- readRDS(path)
        }
      }

      if (length(scsr_results) == 0) {
        stop("No scsrnn objects loaded.")
      }

      # Order names if D3, D7, D14 are present so they look neat
      day_names <- names(scsr_results)
      sort_order <- c("D3", "D7", "D14")
      standard_days <- intersect(sort_order, day_names)
      other_days <- setdiff(day_names, sort_order)
      ordered_names <- c(standard_days, sort_order[sort_order %in% day_names], other_days) |> unique()
      scsr_results <- scsr_results[ordered_names]
      day_names <- names(scsr_results)

      values$scsr_results <- scsr_results

      # Update dynamic UI selectors
      updateSelectInput(session, "explorer_day", choices = day_names, selected = day_names[1])
      updateSelectInput(session, "dist_day1", choices = day_names, selected = day_names[1])
      updateSelectInput(session, "circle_day", choices = day_names, selected = day_names[1])
      updateSelectInput(session, "treemap_day", choices = day_names, selected = day_names[1])
      updateSelectInput(session, "chord_day", choices = day_names, selected = day_names[1])

      if (length(day_names) >= 2) {
        updateSelectInput(session, "dist_day2", choices = day_names, selected = day_names[2])
        updateSelectInput(session, "comp2_day1", choices = day_names, selected = day_names[1])
        updateSelectInput(session, "comp2_day2", choices = day_names, selected = day_names[2])
      } else {
        updateSelectInput(session, "dist_day2", choices = day_names, selected = day_names[1])
      }

      # Update individual visual sliders/inputs to default to load-time values
      updateSliderInput(session, "explorer_lr_threshold", value = input$min_LR_score)
      updateSliderInput(session, "explorer_max_pval", value = input$max_pval)
      updateNumericInput(session, "explorer_min_logFC", value = input$min_logFC)

      updateSliderInput(session, "dist_lr_threshold", value = input$min_LR_score)
      updateSliderInput(session, "dist_max_pval", value = input$max_pval)
      updateNumericInput(session, "dist_min_logFC", value = input$min_logFC)

      updateSliderInput(session, "circle_lr_threshold", value = input$min_LR_score)
      updateSliderInput(session, "circle_max_pval", value = input$max_pval)
      updateNumericInput(session, "circle_min_logFC", value = input$min_logFC)

      updateSliderInput(session, "treemap_lr_threshold", value = input$min_LR_score)
      updateSliderInput(session, "treemap_max_pval", value = input$max_pval)
      updateNumericInput(session, "treemap_min_logFC", value = input$min_logFC)

      updateSliderInput(session, "chord_lr_threshold", value = input$min_LR_score)
      updateSliderInput(session, "chord_max_pval", value = input$max_pval)
      updateNumericInput(session, "chord_min_logFC", value = input$min_logFC)

      updateSliderInput(session, "comp2_lr_threshold", value = input$min_LR_score)
      updateSliderInput(session, "comp2_max_pval", value = input$max_pval)
      updateNumericInput(session, "comp2_min_logFC", value = input$min_logFC)

      updateSliderInput(session, "comp3_lr_threshold", value = input$min_LR_score)
      updateSliderInput(session, "comp3_max_pval", value = input$max_pval)
      updateNumericInput(session, "comp3_min_logFC", value = input$min_logFC)

      values$is_analyzed <- TRUE
      log_message("All datasets loaded and initialized successfully!")
      shiny::showNotification("scsrnn objects loaded successfully!", type = "message")
    }, error = function(e) {
      log_message(paste("CRITICAL LOADING ERROR:", e$message))
      shiny::showNotification(paste("Loading Failed:", e$message), type = "error", duration = NULL)
    }, finally = {
      # 3. Always remove the modal loading overlay to unblock user interaction
      shiny::removeModal()
    })
  })

  # ----------------------------------------------------------------------------
  # D. DYNAMIC DROPDOWN CONTROLS UPDATES
  # ----------------------------------------------------------------------------

  # Tab 2 Explorer selectors
  observe({
    req(values$is_analyzed)
    day <- input$explorer_day
    req(day)
    scsrnn <- values$scsr_results[[day]]
    req(scsrnn)
    clusters <- unique(as.character(slot(scsrnn, "populations")))

    updateSelectInput(session, "explorer_sender", choices = clusters, selected = clusters)
    updateSelectInput(session, "explorer_receiver", choices = clusters, selected = clusters)
  })

  # Tab Chord Diagram selectors
  observe({
    req(values$is_analyzed)
    day <- input$chord_day
    req(day)
    scsrnn <- values$scsr_results[[day]]
    req(scsrnn)
    clusters <- unique(as.character(slot(scsrnn, "populations")))

    updateSelectInput(session, "chord_sender", choices = clusters, selected = clusters)
    updateSelectInput(session, "chord_receiver", choices = clusters, selected = clusters)
  })

  # Select All / None Chord Senders
  observeEvent(input$chord_sender_all, {
    req(values$is_analyzed)
    day <- input$chord_day
    req(day)
    scsrnn <- values$scsr_results[[day]]
    req(scsrnn)
    clusters <- unique(as.character(slot(scsrnn, "populations")))
    updateSelectInput(session, "chord_sender", selected = clusters)
  })

  observeEvent(input$chord_sender_none, {
    updateSelectInput(session, "chord_sender", selected = character(0))
  })

  # Select All / None Chord Receivers
  observeEvent(input$chord_receiver_all, {
    req(values$is_analyzed)
    day <- input$chord_day
    req(day)
    scsrnn <- values$scsr_results[[day]]
    req(scsrnn)
    clusters <- unique(as.character(slot(scsrnn, "populations")))
    updateSelectInput(session, "chord_receiver", selected = clusters)
  })

  observeEvent(input$chord_receiver_none, {
    updateSelectInput(session, "chord_receiver", selected = character(0))
  })

  # Tab 3 Distribution selectors
  observe({
    req(values$is_analyzed)
    day1 <- input$dist_day1
    req(day1)
    scsrnn <- values$scsr_results[[day1]]
    req(scsrnn)
    clusters <- unique(as.character(slot(scsrnn, "populations")))

    updateSelectInput(session, "dist_sender", choices = clusters, selected = clusters[1])
    updateSelectInput(session, "dist_receiver", choices = clusters, selected = if (length(clusters) > 1) clusters[2] else clusters[1])
  })

  # Tab Pathway Circlepack selectors
  observe({
    req(values$is_analyzed)
    day <- input$circle_day
    req(day)
    scsrnn <- values$scsr_results[[day]]
    req(scsrnn)
    clusters <- unique(as.character(slot(scsrnn, "populations")))

    updateSelectInput(session, "circle_sender", choices = clusters, selected = clusters[1])
    updateSelectInput(session, "circle_receiver", choices = clusters, selected = if (length(clusters) > 1) clusters[2] else clusters[1])
  })

  # Tab Pathway Treemap selectors
  observe({
    req(values$is_analyzed)
    day <- input$treemap_day
    req(day)
    scsrnn <- values$scsr_results[[day]]
    req(scsrnn)
    clusters <- unique(as.character(slot(scsrnn, "populations")))

    updateSelectInput(session, "treemap_sender", choices = clusters, selected = clusters[1])
    updateSelectInput(session, "treemap_receiver", choices = clusters, selected = if (length(clusters) > 1) clusters[2] else clusters[1])
  })

  # Tab 4 Compare 2 selectors
  observe({
    req(values$is_analyzed)
    day1 <- input$comp2_day1
    day2 <- input$comp2_day2
    req(day1, day2)

    scsrnn1 <- values$scsr_results[[day1]]
    scsrnn2 <- values$scsr_results[[day2]]
    req(scsrnn1, scsrnn2)
    clusters <- unique(c(
      as.character(slot(scsrnn1, "populations")),
      as.character(slot(scsrnn2, "populations"))
    ))

    updateSelectInput(session, "comp2_sender", choices = clusters, selected = clusters[1])
    updateSelectInput(session, "comp2_receiver", choices = clusters, selected = if (length(clusters) > 1) clusters[2] else clusters[1])
  })

  # Tab 5 Compare 3 selectors
  observe({
    req(values$is_analyzed)
    day_names <- names(values$scsr_results)
    req(length(day_names) >= 3)

    day1 <- day_names[1]
    scsrnn1 <- values$scsr_results[[day1]]
    req(scsrnn1)
    clusters <- unique(as.character(slot(scsrnn1, "populations")))

    updateSelectInput(session, "comp3_sender", choices = clusters, selected = clusters[1])
    updateSelectInput(session, "comp3_receiver", choices = clusters, selected = if (length(clusters) > 1) clusters[2] else clusters[1])
  })

  # Select All / None Senders
  observeEvent(input$sender_all, {
    req(values$is_analyzed)
    day <- input$explorer_day
    req(day)
    scsrnn <- values$scsr_results[[day]]
    req(scsrnn)
    clusters <- unique(as.character(slot(scsrnn, "populations")))
    updateSelectInput(session, "explorer_sender", selected = clusters)
  })

  observeEvent(input$sender_none, {
    updateSelectInput(session, "explorer_sender", selected = character(0))
  })

  # Select All / None Receivers
  observeEvent(input$receiver_all, {
    req(values$is_analyzed)
    day <- input$explorer_day
    req(day)
    scsrnn <- values$scsr_results[[day]]
    req(scsrnn)
    clusters <- unique(as.character(slot(scsrnn, "populations")))
    updateSelectInput(session, "explorer_receiver", selected = clusters)
  })

  observeEvent(input$receiver_none, {
    updateSelectInput(session, "explorer_receiver", selected = character(0))
  })

  # ----------------------------------------------------------------------------
  # E. RENDER TAB 2: GLOBAL LR EXPLORER PLOTS & TABLES
  # ----------------------------------------------------------------------------

  # Filters the global interaction table based on sidebar inputs
  filtered_explorer_table <- reactive({
    req(values$is_analyzed)
    day <- input$explorer_day
    req(day)
    scsrnn <- values$scsr_results[[day]]
    req(scsrnn)

    senders <- input$explorer_sender
    receivers <- input$explorer_receiver

    if (is.null(senders) || length(senders) == 0 || is.null(receivers) || length(receivers) == 0) {
      return(data.frame(
        Sender = character(), Receiver = character(), L = character(), R = character(),
        LR.score = numeric(), L.logFC = numeric(), R.logFC = numeric(), pval = numeric(),
        pathway_name = character(), root_parent = character(), subparent = character(), subsubparent = character()
      ))
    }

    # Fetch safe compiled interactions on the fly!
    df <- get_multiple_interactions_safe(scsrnn, senders, receivers, values$reactome_human)
    if (nrow(df) == 0) {
      return(data.frame(
        Sender = character(), Receiver = character(), L = character(), R = character(),
        LR.score = numeric(), L.logFC = numeric(), R.logFC = numeric(), pval = numeric(),
        pathway_name = character(), root_parent = character(), subparent = character(), subsubparent = character()
      ))
    }

    # Filter by post-hoc LR score
    df <- df |> dplyr::filter(LR.score >= input$explorer_lr_threshold)

    # Filter by P-Value
    df <- df |> dplyr::filter(pval <= input$explorer_max_pval)

    # Filter by Log Fold Change
    df <- df |> dplyr::filter(L.logFC >= input$explorer_min_logFC & R.logFC >= input$explorer_min_logFC)

    # Filter by Gene
    if (nchar(input$explorer_gene) > 0) {
      gene_query <- toupper(str_trim(input$explorer_gene))
      df <- df |> dplyr::filter(
        grepl(gene_query, L) | grepl(gene_query, R)
      )
    }

    return(df)
  })

  # 1. Global Heatmap Reactive Plot Object
  global_heatmap_plot <- reactive({
    req(values$is_analyzed)
    day <- input$explorer_day
    req(day)
    scsrnn <- values$scsr_results[[day]]
    req(scsrnn)

    senders <- input$explorer_sender
    receivers <- input$explorer_receiver

    if (is.null(senders) || length(senders) == 0 || is.null(receivers) || length(receivers) == 0) {
      return(NULL)
    }

    df <- filtered_explorer_table()

    if (is.null(df) || nrow(df) == 0) {
      return(NULL)
    }

    # Count interactions per Sender -> Receiver pair
    df_counts <- df |>
      dplyr::group_by(Sender, Receiver) |>
      dplyr::summarise(Interaction_Count = dplyr::n(), .groups = "drop")

    df_counts <- df_counts |>
      dplyr::mutate(
        Sender = factor(Sender, levels = senders),
        Receiver = factor(Receiver, levels = receivers)
      )

    ggplot(df_counts, aes(x = Sender, y = Receiver, fill = Interaction_Count)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = Interaction_Count), color = "black", fontface = "bold", size = 4.5) +
      scale_fill_gradient(low = "#f8fafc", high = "#6366f1", name = "Interaction Count") +
      labs(
        title = "Intercellular Communication Heatmap",
        subtitle = "Number of ligand-receptor interactions between populations (Sender → Receiver)",
        x = "Sender Populations",
        y = "Receiver Populations"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold", color = "#0f172a"),
        plot.subtitle = element_text(size = 11, color = "#475569"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold", color = "#1e293b"),
        axis.text.y = element_text(size = 11, face = "bold", color = "#1e293b"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
  })

  output$global_heatmap <- renderPlot({
    p <- global_heatmap_plot()
    if (is.null(p)) {
      plot(0, type = "n", axes = FALSE, ann = FALSE)
      text(1, 1, "Please select at least one Sender and Receiver population, and ensure interactions exist.", cex = 1.3, font = 2)
      return(NULL)
    }
    p
  })

  # 1b. Custom Network Bubble Plot Reactive Object
  network_bubble_plot <- reactive({
    req(values$is_analyzed)
    day <- input$explorer_day
    req(day)
    scsrnn <- values$scsr_results[[day]]
    req(scsrnn)

    senders <- input$explorer_sender
    receivers <- input$explorer_receiver

    if (is.null(senders) || length(senders) == 0 || is.null(receivers) || length(receivers) == 0) {
      return(NULL)
    }

    df <- filtered_explorer_table()

    if (is.null(df) || nrow(df) == 0) {
      return(NULL)
    }

    # Count interactions per Sender -> Receiver pair
    df_counts <- df |>
      dplyr::group_by(Sender, Receiver) |>
      dplyr::summarise(Interaction_Count = dplyr::n(), .groups = "drop")

    df_counts <- df_counts |>
      dplyr::mutate(
        Sender = factor(Sender, levels = senders),
        Receiver = factor(Receiver, levels = receivers)
      )

    ggplot(df_counts, aes(x = Sender, y = Receiver, size = Interaction_Count)) +
      geom_point(alpha = 0.85) +
      geom_text(aes(label = Interaction_Count), color = "white", fontface = "bold", size = 4.5) +
      scale_size_continuous(range = c(8, 22), name = "Interaction Count") +
      labs(
        title = "Intercellular Communication Frequencies",
        subtitle = "Bubble size and color represent the number of ligand-receptor interactions (Sender → Receiver)",
        x = "Sender Populations",
        y = "Receiver Populations"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold", color = "#0f172a"),
        plot.subtitle = element_text(size = 11, color = "#475569"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold", color = "#1e293b"),
        axis.text.y = element_text(size = 11, face = "bold", color = "#1e293b"),
        panel.grid.major = element_line(color = "#f1f5f9"),
        panel.grid.minor = element_blank()
      )
  })

  output$network_bubble <- renderPlot({
    p <- network_bubble_plot()
    if (is.null(p)) {
      plot(0, type = "n", axes = FALSE, ann = FALSE)
      text(1, 1, "Please select at least one Sender and Receiver population, and ensure interactions exist.", cex = 1.3, font = 2)
      return(NULL)
    }
    p
  })

  # 2. Bubble Plot View Reactive Object
  global_bubble_plot <- reactive({
    df <- filtered_explorer_table()
    req(df)

    if (nrow(df) == 0) {
      return(NULL)
    }

    # Limit to top 50 interactions to keep the plot clean and readable
    if (nrow(df) > 50) {
      df <- df |>
        dplyr::arrange(dplyr::desc(LR.score)) |>
        dplyr::slice(1:50)
    }

    # Create cleaner labels for plotting
    df <- df |>
      dplyr::mutate(
        Interaction = paste(L, R, sep = " - "),
        Pair = paste(Sender, Receiver, sep = " → ")
      )

    ggplot(df, aes(x = Pair, y = reorder(Interaction, LR.score), size = LR.score)) +
      geom_point(alpha = 0.8) +
      scale_size_continuous(range = c(3.5, 9), name = "LR Score") +
      labs(
        title = "Inferred Directional Cellular Interactions",
        subtitle = "Bubble size and color represent the Regularized LR Score (Sender → Receiver)",
        x = "Cellular Communication Pairs (Sender → Receiver)",
        y = "Ligand - Receptor Interactions"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold", color = "#0f172a"),
        plot.subtitle = element_text(size = 11, color = "#475569"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold", color = "#334155"),
        axis.text.y = element_text(size = 9, color = "#334155"),
        panel.grid.major = element_line(color = "#f1f5f9"),
        panel.grid.minor = element_blank()
      )
  })

  output$global_bubble <- renderPlot({
    df_raw <- filtered_explorer_table()
    if (!is.null(df_raw) && nrow(df_raw) > 50) {
      showNotification("Showing top 50 interactions by LR score to maintain plot readability.", type = "message")
    }
    p <- global_bubble_plot()
    if (is.null(p)) {
      plot(0, type = "n", axes = FALSE, ann = FALSE)
      text(1, 1, "No interactions found for the selected Senders, Receivers, and LR score threshold.", cex = 1.3, font = 2)
      return(NULL)
    }
    p
  })

  # 3. Interactions DT Table
  output$interactions_table <- renderDT({
    df <- filtered_explorer_table()
    req(df)

    # Format for display
    display_df <- df |>
      dplyr::select(Sender, Receiver, L, R, LR.score, L.logFC, R.logFC, pval, pathway_name, root_parent, subparent, subsubparent) |>
      dplyr::mutate(
        LR.score = round(LR.score, 3),
        L.logFC = round(L.logFC, 3),
        R.logFC = round(R.logFC, 3),
        pval = format.pval(pval, digits = 3)
      )

    datatable(display_df,
      filter = "top",
      options = list(pageLength = 10, autoWidth = TRUE, scrollX = TRUE),
      rownames = FALSE
    )
  })

  # Reactive table for Chord Diagram
  filtered_chord_table <- reactive({
    req(values$is_analyzed)
    day <- input$chord_day
    req(day)
    scsrnn <- values$scsr_results[[day]]
    req(scsrnn)

    senders <- input$chord_sender
    receivers <- input$chord_receiver

    if (is.null(senders) || length(senders) == 0 || is.null(receivers) || length(receivers) == 0) {
      return(data.frame())
    }

    # Fetch safe compiled interactions on the fly
    df <- get_multiple_interactions_safe(scsrnn, senders, receivers, values$reactome_human)
    if (nrow(df) == 0) {
      return(data.frame())
    }

    # Filter by user selections
    df <- df |>
      dplyr::filter(LR.score >= input$chord_lr_threshold) |>
      dplyr::filter(pval <= input$chord_max_pval) |>
      dplyr::filter(L.logFC >= input$chord_min_logFC & R.logFC >= input$chord_min_logFC)

    return(df)
  })

  # Render Chord Diagram using circlize
  output$chord_plot <- renderPlot({
    df <- filtered_chord_table()
    req(df)

    if (nrow(df) == 0) {
      return(NULL)
    }

    # Count interactions per Sender -> Receiver pair
    df_counts <- df |>
      dplyr::group_by(Sender, Receiver) |>
      dplyr::summarise(Interaction_Count = dplyr::n(), .groups = "drop")

    if (nrow(df_counts) == 0) {
      return(NULL)
    }

    # Prepare matrix for circlize::chordDiagram
    all_pops <- sort(unique(c(df_counts$Sender, df_counts$Receiver)))
    df_grid <- expand.grid(Sender = all_pops, Receiver = all_pops, stringsAsFactors = FALSE)
    df_complete <- df_grid |>
      dplyr::left_join(df_counts, by = c("Sender", "Receiver")) |>
      dplyr::mutate(Interaction_Count = dplyr::coalesce(Interaction_Count, 0))

    df_wide <- df_complete |>
      tidyr::pivot_wider(names_from = Receiver, values_from = Interaction_Count)

    m <- as.matrix(df_wide[, -1])
    rownames(m) <- df_wide$Sender
    colnames(m) <- df_wide$Sender

    # Assign distinct colors to all unique populations
    groupColors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(length(all_pops))
    names(groupColors) <- all_pops

    # Reset any previous circlize plots to prevent overlapping/clipping
    circos.clear()

    # Render static chordDiagram with directionality and arrows
    circlize::chordDiagram(
      m,
      directional = 1,
      direction.type = c("diffHeight", "arrows"),
      link.arr.type = "big.arrow",
      grid.col = groupColors,
      transparency = 0.4
    )
    circos.clear()
  })

  # ----------------------------------------------------------------------------
  # F. RENDER TAB 3: PATHWAY DISTRIBUTION PLOTS (INTEGRATED COMPARISON)
  # ----------------------------------------------------------------------------

  # Compiles subset of pathways for selected Day 1
  pathway_dist_data1 <- reactive({
    req(values$is_analyzed)
    day <- input$dist_day1
    req(day)

    scsrnn <- values$scsr_results[[day]]
    req(scsrnn)

    df_raw <- get_pathway_frequency_data(
      scsrnn = scsrnn,
      scope = input$dist_scope,
      sender = input$dist_sender,
      receiver = input$dist_receiver,
      reactome_human = values$reactome_human,
      min_lr_score = input$dist_lr_threshold,
      max_pval = input$dist_max_pval,
      min_logFC = input$dist_min_logFC
    )

    if (is.null(df_raw) || nrow(df_raw) == 0 || is.null(values$ref_couleurs)) {
      return(NULL)
    }

    # Group and count occurrences
    concat_csv_subparent_root <- df_raw |>
      dplyr::group_by(root_parent, subparent) |>
      dplyr::count() |>
      dplyr::ungroup()

    if (nrow(concat_csv_subparent_root) == 0) {
      return(NULL)
    }

    # Join with our fixed global colors lookup
    df_couleur <- concat_csv_subparent_root |>
      dplyr::left_join(values$ref_couleurs, by = c("root_parent", "subparent"))

    return(df_couleur)
  })

  # Compiles subset of pathways for selected Day 2
  pathway_dist_data2 <- reactive({
    req(values$is_analyzed)
    day <- input$dist_day2
    req(day)

    scsrnn <- values$scsr_results[[day]]
    req(scsrnn)

    df_raw <- get_pathway_frequency_data(
      scsrnn = scsrnn,
      scope = input$dist_scope,
      sender = input$dist_sender,
      receiver = input$dist_receiver,
      reactome_human = values$reactome_human,
      min_lr_score = input$dist_lr_threshold,
      max_pval = input$dist_max_pval,
      min_logFC = input$dist_min_logFC
    )

    if (is.null(df_raw) || nrow(df_raw) == 0 || is.null(values$ref_couleurs)) {
      return(NULL)
    }

    # Group and count occurrences
    concat_csv_subparent_root <- df_raw |>
      dplyr::group_by(root_parent, subparent) |>
      dplyr::count() |>
      dplyr::ungroup()

    if (nrow(concat_csv_subparent_root) == 0) {
      return(NULL)
    }

    # Join with our fixed global colors lookup
    df_couleur <- concat_csv_subparent_root |>
      dplyr::left_join(values$ref_couleurs, by = c("root_parent", "subparent"))

    return(df_couleur)
  })

  # Compiles subset of pathways for Circlepack
  pathway_circle_data <- reactive({
    req(values$is_analyzed)
    day <- input$circle_day
    req(day)

    scsrnn <- values$scsr_results[[day]]
    req(scsrnn)

    df_raw <- get_pathway_frequency_data(
      scsrnn = scsrnn,
      scope = input$circle_scope,
      sender = input$circle_sender,
      receiver = input$circle_receiver,
      reactome_human = values$reactome_human,
      min_lr_score = input$circle_lr_threshold,
      max_pval = input$circle_max_pval,
      min_logFC = input$circle_min_logFC
    )

    if (is.null(df_raw) || nrow(df_raw) == 0 || is.null(values$ref_couleurs)) {
      return(NULL)
    }

    # Group and count occurrences
    concat_csv_subparent_root <- df_raw |>
      dplyr::group_by(root_parent, subparent) |>
      dplyr::count() |>
      dplyr::ungroup()

    return(concat_csv_subparent_root)
  })

  # Compiles subset of pathways for Treemap
  pathway_treemap_data <- reactive({
    req(values$is_analyzed)
    day <- input$treemap_day
    req(day)

    scsrnn <- values$scsr_results[[day]]
    req(scsrnn)

    df_raw <- get_pathway_frequency_data(
      scsrnn = scsrnn,
      scope = input$treemap_scope,
      sender = input$treemap_sender,
      receiver = input$treemap_receiver,
      reactome_human = values$reactome_human,
      min_lr_score = input$treemap_lr_threshold,
      max_pval = input$treemap_max_pval,
      min_logFC = input$treemap_min_logFC
    )

    if (is.null(df_raw) || nrow(df_raw) == 0 || is.null(values$ref_couleurs)) {
      return(NULL)
    }

    # Group and count occurrences
    concat_csv_subparent_root <- df_raw |>
      dplyr::group_by(root_parent, subparent) |>
      dplyr::count() |>
      dplyr::ungroup()

    return(concat_csv_subparent_root)
  })

  # Dynamic Headers
  output$header_day1 <- renderUI({
    req(input$dist_day1)
    span(icon("chart-bar"), paste(input$dist_day1, "Pathway Frequencies"))
  })

  output$header_day2 <- renderUI({
    req(input$dist_day2)
    span(icon("chart-bar"), paste(input$dist_day2, "Pathway Frequencies"))
  })

  # 1a. Stacked Subparents Bar Plot - Day 1 (Left) - Reactive ggplot Object
  pathway_dist_plot1_ggplot <- reactive({
    df <- pathway_dist_data1()
    if (is.null(df) || nrow(df) == 0) {
      return(NULL)
    }

    if (input$dist_level == "root_parent") {
      # Group by root parent only
      df_root <- df |>
        dplyr::group_by(root_parent) |>
        dplyr::summarise(n = sum(n), .groups = "drop")

      p <- ggplot(df_root, aes(
        x = reorder(root_parent, n), y = n,
        fill = root_parent,
        text = paste0(
          "<b>Root Parent :</b> ", root_parent, "<br>",
          "<b>Total (n) :</b> ", n
        )
      )) +
        geom_col() +
        coord_flip() +
        scale_fill_viridis_d(option = "mako", guide = "none") +
        theme_minimal() +
        theme(
          legend.position = "none",
          panel.grid.minor = element_blank(),
          axis.text.y = element_text(size = 9)
        ) +
        labs(
          x = "Root Parent Pathway",
          y = "Interaction Frequencies (n)"
        )
    } else {
      # Subparent level (matches root_parents_distribution_per_day.R)
      p <- ggplot(df, aes(
        x = root_parent, y = n,
        fill = as.factor(index_couleur),
        text = paste0(
          "<b>Subparent :</b> ", subparent, "<br>",
          "<b>Total (n) :</b> ", n
        )
      )) +
        geom_col(position = position_stack(reverse = FALSE)) +
        coord_flip() +
        scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(values$max_colors)) +
        theme_minimal() +
        theme(
          legend.position = "none",
          panel.grid.minor = element_blank(),
          axis.text.y = element_text(size = 9)
        ) +
        labs(
          x = "Root Parent Pathway",
          y = "Interaction Frequencies (n)"
        )
    }
    p
  })

  # 1a. Stacked Subparents Bar Plot (Plotly) - Day 1 (Left)
  output$pathway_stacked_bar1 <- renderPlotly({
    p <- pathway_dist_plot1_ggplot()
    if (is.null(p)) {
      return(plotly_empty(type = "scatter", mode = "markers") |>
        layout(title = "No pathway mappings found."))
    }
    ggplotly(p, tooltip = c("text")) |>
      layout(showlegend = FALSE)
  })

  # 1b. Stacked Subparents Bar Plot - Day 2 (Right) - Reactive ggplot Object
  pathway_dist_plot2_ggplot <- reactive({
    df <- pathway_dist_data2()
    if (is.null(df) || nrow(df) == 0) {
      return(NULL)
    }

    if (input$dist_level == "root_parent") {
      # Group by root parent only
      df_root <- df |>
        dplyr::group_by(root_parent) |>
        dplyr::summarise(n = sum(n), .groups = "drop")

      p <- ggplot(df_root, aes(
        x = reorder(root_parent, n), y = n,
        fill = root_parent,
        text = paste0(
          "<b>Root Parent :</b> ", root_parent, "<br>",
          "<b>Total (n) :</b> ", n
        )
      )) +
        geom_col() +
        coord_flip() +
        scale_fill_viridis_d(option = "mako", guide = "none") +
        theme_minimal() +
        theme(
          legend.position = "none",
          panel.grid.minor = element_blank(),
          axis.text.y = element_text(size = 9)
        ) +
        labs(
          x = "Root Parent Pathway",
          y = "Interaction Frequencies (n)"
        )
    } else {
      # Subparent level (matches root_parents_distribution_per_day.R)
      p <- ggplot(df, aes(
        x = root_parent, y = n,
        fill = as.factor(index_couleur),
        text = paste0(
          "<b>Subparent :</b> ", subparent, "<br>",
          "<b>Total (n) :</b> ", n
        )
      )) +
        geom_col(position = position_stack(reverse = FALSE)) +
        coord_flip() +
        scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(values$max_colors)) +
        theme_minimal() +
        theme(
          legend.position = "none",
          panel.grid.minor = element_blank(),
          axis.text.y = element_text(size = 9)
        ) +
        labs(
          x = "Root Parent Pathway",
          y = "Interaction Frequencies (n)"
        )
    }
    p
  })

  # 1b. Stacked Subparents Bar Plot (Plotly) - Day 2 (Right)
  output$pathway_stacked_bar2 <- renderPlotly({
    p <- pathway_dist_plot2_ggplot()
    if (is.null(p)) {
      return(plotly_empty(type = "scatter", mode = "markers") |>
        layout(title = "No pathway mappings found."))
    }
    ggplotly(p, tooltip = c("text")) |>
      layout(showlegend = FALSE)
  })

  # 2. Circlepack Hierarchy Graph (circlepackeR) - dedicated to Circlepack
  output$pathway_circlepack <- renderCirclepackeR({
    df <- pathway_circle_data()
    if (is.null(df) || nrow(df) == 0) {
      return(NULL)
    }

    # Build Edges: Racine -> Root Parent -> Subparent
    edges <- data.frame(
      from = c(rep("Racine", times = length(unique(df$root_parent))), df$root_parent),
      to = c(unique(df$root_parent), df$subparent)
    ) |> distinct()

    vertices <- data.frame(
      name = c("Racine", unique(df$root_parent), unique(df$subparent))
    ) |> distinct()

    # Map sizes/weights
    subparent_weights <- df |>
      group_by(subparent) |>
      summarise(valeur = sum(n), .groups = "drop")

    vertices <- vertices |>
      left_join(subparent_weights, by = c("name" = "subparent")) |>
      mutate(valeur = ifelse(is.na(valeur), 0, valeur))

    # dataset
    edges_no_root <- edges |>
      filter(from != "Racine")
    group <- edges_no_root$from
    subgroup <- edges_no_root$to
    value <- edges_no_root |>
      left_join(vertices, by = c("to" = "name")) |>
      pull(valeur)

    data <- data.frame(group, subgroup, value)
    data$pathString <- paste("world", data$group, data$subgroup, sep = "/")
    population <- data.tree::as.Node(data)

    p <- circlepackeR::circlepackeR(population, size = "value", color_min = "hsl(56,80%,80%)", color_max = "hsl(341,30%,40%)")
    p <- htmlwidgets::onRender(p, "
    function(el, x) {
        d3.select(el).selectAll('circle').filter(function(d) { return d.depth === 0; }).style('display', 'none');
    }
    ")
    return(p)
  })

  output$pathway_treemap <- renderD3tree2({
    df <- pathway_treemap_data()
    if (is.null(df) || nrow(df) == 0) {
      return(NULL)
    }

    # Group and count (pathway_treemap_data returns root_parent, subparent, n)
    data <- df |>
      dplyr::rename(group = root_parent, subgroup = subparent, value = n) |>
      dplyr::select(group, subgroup, value) |>
      dplyr::distinct()

    if (nrow(data) == 0) {
      return(NULL)
    }

    # Map colors consistently based on root parent names
    unique_groups <- sort(unique(data$group))
    unique_groups <- unique_groups[unique_groups != "" & !is.na(unique_groups)]

    if (!is.null(values$root_colors)) {
      mapped_colors <- values$root_colors[unique_groups]
    } else {
      # Fallback dynamic palette if database mapping is not loaded yet
      mapped_colors <- colorRampPalette(c(
        RColorBrewer::brewer.pal(12, "Set3"),
        RColorBrewer::brewer.pal(8, "Set2")
      ))(length(unique_groups))
      names(mapped_colors) <- unique_groups
    }

    pdf(NULL)

    p <- treemap::treemap(data,
      index = c("group", "subgroup"),
      vSize = "value",
      type = "index",
      palette = mapped_colors,
      bg.labels = c("white"),
      align.labels = list(c("center", "center"), c("right", "bottom")),
      draw = FALSE
    )

    dev.off()

    inter <- d3treeR::d3tree2(p, rootname = "Main pathways")

    return(inter)
  })

  # ----------------------------------------------------------------------------
  # G. RENDER TAB 4: COMPARE UNIQUE INTERACTIONS (2 DAYS)
  # ----------------------------------------------------------------------------
  # Reactive ggplot Object for Compare Unique Interactions (2 Days)
  compare2_ggplot <- reactive({
    req(values$is_analyzed)
    day1 <- input$comp2_day1
    day2 <- input$comp2_day2
    sender <- input$comp2_sender
    receiver <- input$comp2_receiver
    column_name <- input$comp2_level
    req(day1, day2, sender, receiver)

    if (day1 == day2) {
      return(NULL)
    }

    # 1. Fetch safe pair interactions dynamically for day 1 and day 2
    scsrnn1 <- values$scsr_results[[day1]]
    scsrnn2 <- values$scsr_results[[day2]]
    req(scsrnn1, scsrnn2)

    df_day_1_raw <- get_pair_interactions_safe(scsrnn1, sender, receiver, values$reactome_human)
    if (nrow(df_day_1_raw) > 0) {
      df_day_1_raw <- df_day_1_raw |>
        dplyr::filter(LR.score >= input$comp2_lr_threshold) |>
        dplyr::filter(pval <= input$comp2_max_pval) |>
        dplyr::filter(L.logFC >= input$comp2_min_logFC & R.logFC >= input$comp2_min_logFC)
    }

    df_day_2_raw <- get_pair_interactions_safe(scsrnn2, sender, receiver, values$reactome_human)
    if (nrow(df_day_2_raw) > 0) {
      df_day_2_raw <- df_day_2_raw |>
        dplyr::filter(LR.score >= input$comp2_lr_threshold) |>
        dplyr::filter(pval <= input$comp2_max_pval) |>
        dplyr::filter(L.logFC >= input$comp2_min_logFC & R.logFC >= input$comp2_min_logFC)
    }

    # Process Day 1 dataframe
    if (nrow(df_day_1_raw) > 0) {
      df_day_1 <- df_day_1_raw |>
        mutate(across(where(is.character), ~ if_else(.x == "", NA_character_, .x))) |>
        drop_na(pathway_name) |>
        unite(col = "unique_id", L, R, sep = "_", remove = FALSE)
    } else {
      df_day_1 <- data.frame(unique_id = character())
    }

    # Process Day 2 dataframe
    if (nrow(df_day_2_raw) > 0) {
      df_day_2 <- df_day_2_raw |>
        mutate(across(where(is.character), ~ if_else(.x == "", NA_character_, .x))) |>
        drop_na(pathway_name) |>
        unite(col = "unique_id", L, R, sep = "_", remove = FALSE)
    } else {
      df_day_2 <- data.frame(unique_id = character())
    }

    if (nrow(df_day_1) == 0 && nrow(df_day_2) == 0) {
      return(NULL)
    }

    # Find non-common interactions
    non_common_interactions <- append(
      setdiff(df_day_1$unique_id, df_day_2$unique_id),
      setdiff(df_day_2$unique_id, df_day_1$unique_id)
    )

    df_non_common_day1 <- df_day_1 |> dplyr::filter(unique_id %in% non_common_interactions)
    df_non_common_day2 <- df_day_2 |> dplyr::filter(unique_id %in% non_common_interactions)

    if (nrow(df_non_common_day1) == 0 && nrow(df_non_common_day2) == 0) {
      return(NULL)
    }

    # Gather counts for selected hierarchy column
    df_day_1_pathways <- data.frame()
    if (nrow(df_non_common_day1) > 0 && column_name %in% colnames(df_non_common_day1)) {
      all_pathways_day_1 <- unlist(strsplit(df_non_common_day1[[column_name]], "; "))
      if (length(all_pathways_day_1) > 0) {
        df_day_1_pathways <- data.frame(pathways = all_pathways_day_1) |>
          count(pathways, sort = TRUE) |>
          rename(count = n) |>
          mutate(day = day1)
      }
    }

    df_day_2_pathways <- data.frame()
    if (nrow(df_non_common_day2) > 0 && column_name %in% colnames(df_non_common_day2)) {
      all_pathways_day_2 <- unlist(strsplit(df_non_common_day2[[column_name]], "; "))
      if (length(all_pathways_day_2) > 0) {
        df_day_2_pathways <- data.frame(pathways = all_pathways_day_2) |>
          count(pathways, sort = TRUE) |>
          rename(count = n) |>
          mutate(day = day2)
      }
    }

    if (nrow(df_day_1_pathways) == 0 && nrow(df_day_2_pathways) == 0) {
      return(NULL)
    }

    # Combine lists and format for diverging columns
    df_combined <- bind_rows(df_day_1_pathways, df_day_2_pathways) |>
      pivot_wider(names_from = day, values_from = count, values_fill = 0)

    # Insure both day columns exist
    if (!day1 %in% colnames(df_combined)) df_combined[[day1]] <- 0
    if (!day2 %in% colnames(df_combined)) df_combined[[day2]] <- 0

    df_combined <- df_combined |>
      mutate(Sort_Score = .data[[day1]] - .data[[day2]]) |>
      arrange(Sort_Score) |>
      mutate(pathways = factor(pathways, levels = pathways)) |>
      pivot_longer(cols = c(all_of(day1), all_of(day2)), names_to = "day", values_to = "Count") |>
      mutate(Plot_Count = ifelse(day == day1, -Count, Count))

    max_abs_count <- max(abs(df_combined$Plot_Count), 1)

    # Render ggplot diverging chart
    ggplot(df_combined, aes(x = pathways, y = Plot_Count, fill = day)) +
      geom_col(width = 0.75, color = "white", linewidth = 0.2) +
      coord_flip() +
      scale_fill_manual(values = setNames(c("#4682B4", "#E15759"), c(day1, day2))) +
      scale_y_continuous(labels = abs, breaks = seq(-max_abs_count, max_abs_count, by = max(1, round(max_abs_count / 10)))) +
      labs(
        title = glue("Pathway Comparison from {sender} to {receiver} ({day1} vs {day2})"),
        subtitle = glue("Hierarchy: {column_name} (Shows counts of unique, non-common interactions)"),
        y = "Interaction Counts",
        x = "Pathways"
      ) +
      theme_minimal() +
      theme(
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 9, color = "#334155")
      )
  })

  output$compare2_plot <- renderPlot({
    p <- compare2_ggplot()
    if (is.null(p)) {
      plot(0, type = "n", axes = FALSE, ann = FALSE)
      if (input$comp2_day1 == input$comp2_day2) {
        text(1, 1, "Please select two distinct datasets to compare.", cex = 1.5)
      } else {
        text(1, 1, "No unique/non-common interactions found or no mapped pathways.", cex = 1.5)
      }
      return(NULL)
    }
    p
  })

  # ----------------------------------------------------------------------------
  # H. RENDER TAB 5: COMPARE 3 DATASETS (DUMBBELL DOT PLOT)
  # ----------------------------------------------------------------------------
  compare3_ggplot <- reactive({
    req(values$is_analyzed)
    day_names <- names(values$scsr_results)

    # App requires exactly or at least 3 datasets for the dumbbell comparison
    if (length(day_names) < 3) {
      return(NULL)
    }

    day_1 <- day_names[1]
    day_2 <- day_names[2]
    day_3 <- day_names[3]

    sender <- input$comp3_sender
    receiver <- input$comp3_receiver
    column_name <- input$comp3_level
    req(sender, receiver, column_name)

    scsrnn1 <- values$scsr_results[[day_1]]
    scsrnn2 <- values$scsr_results[[day_2]]
    scsrnn3 <- values$scsr_results[[day_3]]
    req(scsrnn1, scsrnn2, scsrnn3)

    df1_raw <- get_pair_interactions_safe(scsrnn1, sender, receiver, values$reactome_human)
    if (nrow(df1_raw) > 0) {
      df1_raw <- df1_raw |>
        dplyr::filter(LR.score >= input$comp3_lr_threshold) |>
        dplyr::filter(pval <= input$comp3_max_pval) |>
        dplyr::filter(L.logFC >= input$comp3_min_logFC & R.logFC >= input$comp3_min_logFC)
    }

    df2_raw <- get_pair_interactions_safe(scsrnn2, sender, receiver, values$reactome_human)
    if (nrow(df2_raw) > 0) {
      df2_raw <- df2_raw |>
        dplyr::filter(LR.score >= input$comp3_lr_threshold) |>
        dplyr::filter(pval <= input$comp3_max_pval) |>
        dplyr::filter(L.logFC >= input$comp3_min_logFC & R.logFC >= input$comp3_min_logFC)
    }

    df3_raw <- get_pair_interactions_safe(scsrnn3, sender, receiver, values$reactome_human)
    if (nrow(df3_raw) > 0) {
      df3_raw <- df3_raw |>
        dplyr::filter(LR.score >= input$comp3_lr_threshold) |>
        dplyr::filter(pval <= input$comp3_max_pval) |>
        dplyr::filter(L.logFC >= input$comp3_min_logFC & R.logFC >= input$comp3_min_logFC)
    }

    # 1. Process data for Day 1
    if (nrow(df1_raw) > 0) {
      df_day_1 <- df1_raw |>
        mutate(across(where(is.character), ~ if_else(.x == "", NA_character_, .x))) |>
        drop_na(pathway_name) |>
        unite(col = "unique_id", L, R, sep = "_", remove = FALSE)
    } else {
      df_day_1 <- data.frame(unique_id = character())
    }

    # 2. Process data for Day 2
    if (nrow(df2_raw) > 0) {
      df_day_2 <- df2_raw |>
        mutate(across(where(is.character), ~ if_else(.x == "", NA_character_, .x))) |>
        drop_na(pathway_name) |>
        unite(col = "unique_id", L, R, sep = "_", remove = FALSE)
    } else {
      df_day_2 <- data.frame(unique_id = character())
    }

    # 3. Process data for Day 3
    if (nrow(df3_raw) > 0) {
      df_day_3 <- df3_raw |>
        mutate(across(where(is.character), ~ if_else(.x == "", NA_character_, .x))) |>
        drop_na(pathway_name) |>
        unite(col = "unique_id", L, R, sep = "_", remove = FALSE)
    } else {
      df_day_3 <- data.frame(unique_id = character())
    }

    # Find three-way intersect
    common_rows <- intersect(df_day_1$unique_id, intersect(df_day_2$unique_id, df_day_3$unique_id))

    # Subtract common rows
    df_day_1 <- df_day_1 |> dplyr::filter(!unique_id %in% common_rows)
    df_day_2 <- df_day_2 |> dplyr::filter(!unique_id %in% common_rows)
    df_day_3 <- df_day_3 |> dplyr::filter(!unique_id %in% common_rows)

    # Create counts
    get_pathway_counts <- function(df, day_lbl) {
      if (nrow(df) == 0 || !column_name %in% colnames(df)) {
        return(data.frame())
      }
      all_paths <- df |>
        pull(.data[[column_name]]) |>
        strsplit(split = "; ", fixed = TRUE) |>
        unlist()
      if (length(all_paths) == 0) {
        return(data.frame())
      }

      data.frame(pathways = all_paths) |>
        count(pathways, sort = TRUE) |>
        rename(count = n) |>
        mutate(day = day_lbl)
    }

    pathways_day_1 <- get_pathway_counts(df_day_1, day_1)
    pathways_day_2 <- get_pathway_counts(df_day_2, day_2)
    pathways_day_3 <- get_pathway_counts(df_day_3, day_3)

    # Bind rows dynamically based on data presence
    df_list <- list(pathways_day_1, pathways_day_2, pathways_day_3)
    df_list <- df_list[sapply(df_list, nrow) > 0]

    if (length(df_list) == 0) {
      return(NULL)
    }

    df_combined_long <- bind_rows(df_list)
    days_in_dataset <- unique(df_combined_long$day)

    # Sorting order based on spreads
    sorting_order <- df_combined_long |>
      arrange(day) |>
      pivot_wider(names_from = day, values_from = count, values_fill = 0)

    # Make sure all days exist in pivoting
    for (d_name in days_in_dataset) {
      if (!d_name %in% colnames(sorting_order)) sorting_order[[d_name]] <- 0
    }

    sorting_order <- sorting_order |>
      mutate(Spread = do.call(pmax, c(across(all_of(days_in_dataset)), na.rm = TRUE)) -
        do.call(pmin, c(across(all_of(days_in_dataset)), na.rm = TRUE))) |>
      arrange(Spread) |>
      pull(pathways)

    df_plot_dumbbell <- df_combined_long |>
      mutate(
        pathways = factor(pathways, levels = sorting_order),
        day = factor(day, levels = days_in_dataset)
      )

    # Adjust overlaps slightly for beautiful rendering (v_offset)
    df_plot_adjusted <- df_plot_dumbbell |>
      group_by(pathways, count) |>
      mutate(
        v_offset = if (n() > 1) {
          seq(-0.14, 0.14, length.out = n())
        } else {
          0
        }
      ) |>
      ungroup() |>
      mutate(pathway_numeric = as.numeric(as.factor(pathways)) + v_offset)

    # Draw ggplot dumbbell plot
    ggplot(df_plot_adjusted) +
      geom_point(aes(x = count, y = pathway_numeric, color = day),
        size = 4, alpha = 0.85
      ) +
      scale_color_manual(values = setNames(c("#2CA02C", "#4682B4", "#E15759"), c(day_1, day_2, day_3))) +
      labs(
        title = glue("Pathway Comparison from {sender} to {receiver}"),
        subtitle = glue("Days: {day_1}, {day_2}, {day_3} | Hierarchy: {column_name} (Intersection excluded)"),
        x = "Interaction Frequency Count",
        y = "Reactome Pathways",
        color = "Condition / Day"
      ) +
      theme_minimal() +
      scale_x_continuous(breaks = seq(0, max(df_plot_adjusted$count) + 1, by = 1)) +
      scale_y_continuous(
        breaks = seq_along(levels(df_plot_adjusted$pathways)),
        labels = levels(df_plot_adjusted$pathways)
      ) +
      theme(
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 9, color = "#334155")
      )
  })

  output$compare3_plot <- renderPlot({
    p <- compare3_ggplot()
    if (is.null(p)) {
      plot(0, type = "n", axes = FALSE, ann = FALSE)
      day_names <- names(values$scsr_results)
      if (length(day_names) < 3) {
        text(1, 1, "Tab requires at least 3 datasets. Select a Splitting column to create subsets.", cex = 1.5)
      } else {
        text(1, 1, "No unique pathway counts available across all days.", cex = 1.5)
      }
      return(NULL)
    }
    p
  })

  # --- PNG DOWNLOAD HANDLERS FOR ALL PLOTS ---

  output$download_heatmap_png <- downloadHandler(
    filename = function() { paste0("global_heatmap-", Sys.Date(), ".png") },
    content = function(file) {
      p <- global_heatmap_plot()
      if (!is.null(p)) ggsave(file, plot = p, device = "png", width = 10, height = 8, dpi = 300)
    }
  )

  output$download_network_bubble_png <- downloadHandler(
    filename = function() { paste0("network_bubble-", Sys.Date(), ".png") },
    content = function(file) {
      p <- network_bubble_plot()
      if (!is.null(p)) ggsave(file, plot = p, device = "png", width = 10, height = 8, dpi = 300)
    }
  )

  output$download_global_bubble_png <- downloadHandler(
    filename = function() { paste0("global_bubble-", Sys.Date(), ".png") },
    content = function(file) {
      p <- global_bubble_plot()
      if (!is.null(p)) ggsave(file, plot = p, device = "png", width = 11, height = 9, dpi = 300)
    }
  )

  output$download_chord_png <- downloadHandler(
    filename = function() { paste0("chord_diagram-", Sys.Date(), ".png") },
    content = function(file) {
      png(file, width = 10, height = 10, units = "in", res = 300)
      df <- filtered_chord_table()
      if (!is.null(df) && nrow(df) > 0) {
        df_counts <- df |>
          dplyr::group_by(Sender, Receiver) |>
          dplyr::summarise(Interaction_Count = dplyr::n(), .groups = "drop")
        
        if (nrow(df_counts) > 0) {
          all_pops <- sort(unique(c(df_counts$Sender, df_counts$Receiver)))
          df_grid <- expand.grid(Sender = all_pops, Receiver = all_pops, stringsAsFactors = FALSE)
          df_complete <- df_grid |>
            dplyr::left_join(df_counts, by = c("Sender", "Receiver")) |>
            dplyr::mutate(Interaction_Count = dplyr::coalesce(Interaction_Count, 0))

          df_wide <- df_complete |>
            tidyr::pivot_wider(names_from = Receiver, values_from = Interaction_Count)

          m <- as.matrix(df_wide[, -1])
          rownames(m) <- df_wide$Sender
          colnames(m) <- df_wide$Sender

          groupColors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(length(all_pops))
          names(groupColors) <- all_pops

          circos.clear()
          circlize::chordDiagram(
            m,
            directional = 1,
            direction.type = c("diffHeight", "arrows"),
            link.arr.type = "big.arrow",
            grid.col = groupColors,
            transparency = 0.4
          )
          circos.clear()
        }
      }
      dev.off()
    }
  )

  output$download_stacked1_png <- downloadHandler(
    filename = function() { paste0("pathway_stacked_bar1-", Sys.Date(), ".png") },
    content = function(file) {
      p <- pathway_dist_plot1_ggplot()
      if (!is.null(p)) ggsave(file, plot = p, device = "png", width = 10, height = 8, dpi = 300)
    }
  )

  output$download_stacked2_png <- downloadHandler(
    filename = function() { paste0("pathway_stacked_bar2-", Sys.Date(), ".png") },
    content = function(file) {
      p <- pathway_dist_plot2_ggplot()
      if (!is.null(p)) ggsave(file, plot = p, device = "png", width = 10, height = 8, dpi = 300)
    }
  )

  output$download_compare2_png <- downloadHandler(
    filename = function() { paste0("compare2_plot-", Sys.Date(), ".png") },
    content = function(file) {
      p <- compare2_ggplot()
      if (!is.null(p)) ggsave(file, plot = p, device = "png", width = 10, height = 8, dpi = 300)
    }
  )

  output$download_compare3_png <- downloadHandler(
    filename = function() { paste0("compare3_plot-", Sys.Date(), ".png") },
    content = function(file) {
      p <- compare3_ggplot()
      if (!is.null(p)) ggsave(file, plot = p, device = "png", width = 10, height = 9, dpi = 300)
    }
  )
}

# ==============================================================================
# 5. RUN SHINY APPLICATION
# ==============================================================================
shinyApp(ui = ui, server = server)
