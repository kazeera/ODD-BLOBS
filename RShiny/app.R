# app.R ---- ODD-BLOBS Shiny visualizer (flex channel mapping, autorun, PDF export) ----

library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(grDevices)
library(bslib)
library(viridisLite)

# ---- Constants expected by functions4.R ----
OFF <- 0L
ON  <- 1L
FORK_OPEN  <- 2L
FORK_CLOSE <- 3L

# ---- Source ODD-BLOBS functions ----
source("functions4.R")

# ---- Robust-ish reader (handles UTF-16LE, incomplete final line) ----
read_fiber_table <- function(path) {
  raw4 <- readBin(path, "raw", n = 4)
  is_utf16le <- length(raw4) >= 2 && raw4[1] == as.raw(0xFF) && raw4[2] == as.raw(0xFE)
  enc <- if (is_utf16le) "UTF-16LE" else "UTF-8"
  
  df <- suppressWarnings(
    tryCatch(
      read.table(
        path,
        sep = "\t",
        header = TRUE,
        fileEncoding = enc,
        fill = TRUE,
        comment.char = "",
        check.names = FALSE
      ),
      error = function(e) NULL
    )
  )
  
  if (is.null(df) || ncol(df) < 2) {
    df <- suppressWarnings(
      read.table(
        path,
        sep = "\t",
        header = TRUE,
        fill = TRUE,
        comment.char = "",
        check.names = FALSE
      )
    )
  }
  
  df
}

# ---- Fix: some versions of get_prot_percents count FORK_OPEN twice ----
get_prot_percents_fixed <- function(prot_array, tract_array) {
  prot_in_tracts <- count_matches(prot_array, tract_array, ON)
  prot_in_forks  <- count_matches(prot_array, tract_array, FORK_OPEN) +
    count_matches(prot_array, tract_array, FORK_CLOSE)
  prot_in_unrep  <- count_matches(prot_array, tract_array, OFF)
  total_prot_ONS <- sum(prot_array == ON)
  
  prot_counts <- c(prot_in_tracts, prot_in_forks, prot_in_unrep, total_prot_ONS)
  if (prot_counts[4] == 0) return(c(0, 0, 0, 0))
  prot_counts * 100 / prot_counts[4]
}


# ---- Helper: normalize hex color strings ----
normalize_hex <- function(x, default) {
  if (is.null(x) || is.na(x)) return(default)
  x <- trimws(x)
  if (grepl("^#[0-9A-Fa-f]{6}$", x)) return(x)
  if (grepl("^[0-9A-Fa-f]{6}$", x)) return(paste0("#", x))
  default
}

# ---- Runner: uses user-chosen mapping for DNA/BrdU/Prot1/Prot2 ----
run_oddblobs_with_mapping <- function(df,
                                      dna_idx, brdu_idx, p1_idx, p2_idx,
                                      brdu_thres, p1_thres, p2_thres,
                                      pR, pU, smooth_val) {
  
  channel_cols <- which(grepl("^Channel", names(df)))
  if (length(channel_cols) < 3) stop("Need at least 3 columns named 'Channel ...'.")
  
  get_chan <- function(idx) {
    if (is.na(idx) || idx < 1 || idx > length(channel_cols)) return(NULL)
    suppressWarnings(as.numeric(df[[channel_cols[idx]]]))
  }
  
  dna  <- get_chan(dna_idx)
  brdu <- get_chan(brdu_idx)
  p1   <- get_chan(p1_idx)
  p2   <- get_chan(p2_idx)
  
  if (is.null(dna) || is.null(brdu) || is.null(p1)) stop("DNA, BrdU, and Protein 1 must be selected.")
  
  dna[is.na(dna)] <- 0
  brdu[is.na(brdu)] <- 0
  p1[is.na(p1)] <- 0
  has_p2 <- !is.null(p2)
  if (has_p2) p2[is.na(p2)] <- 0
  
  brdu_thres <- as.integer(brdu_thres)
  p1_thres <- as.integer(p1_thres)
  p2_thres <- as.integer(p2_thres)
  pR <- as.integer(pR)
  pU <- as.integer(pU)
  smooth_val <- as.integer(smooth_val)
  
  # BrdU -> tracts/forks
  tract_array1 <- threshold_array(brdu, brdu_thres)
  tract_array2 <- smooth_it(tract_array1, smooth_val)
  
  tracts_df <- find_tracts(tract_array2)
  tract_starts  <- tracts_df[["Starts.At"]]
  tract_ends    <- tracts_df[["Ends.At"]]
  tract_lengths <- tracts_df[["Length"]]
  
  tract_array3 <- add_forks(
    tract_array2,
    tract_starts, tract_ends, tract_lengths,
    pU, pR,
    length(tract_array2)
  )
  
  table1 <- find_regions(tract_array3)
  
  # Proteins
  p1_bin <- smooth_it(threshold_array(p1, p1_thres), smooth_val)
  table1 <- cbind(
    table1,
    Protein1 = prot_indices_in_region(p1_bin, nrow(table1), table1$Start, table1$End)
  )
  
  table2 <- get_region_sizes(tract_array3)
  table2 <- cbind(
    table2,
    Prot1Percent = get_prot_percents_fixed(p1_bin, tract_array3)
  )
  
  if (has_p2) {
    p2_bin <- smooth_it(threshold_array(p2, p2_thres), smooth_val)
    
    table1 <- cbind(
      table1,
      Protein2 = prot_indices_in_region(p2_bin, nrow(table1), table1$Start, table1$End)
    )
    
    coloc <- as.integer(p1_bin == ON & p2_bin == ON)
    table1 <- cbind(
      table1,
      Colocalization = prot_indices_in_region(coloc, nrow(table1), table1$Start, table1$End)
    )
    
    table2 <- cbind(
      table2,
      Prot2Percent = get_prot_percents_fixed(p2_bin, tract_array3),
      ColocalizationPercent = get_prot_percents_fixed(coloc, tract_array3)
    )
  }
  
  table1$Code <- NULL
  table2 <- table2 %>%
    mutate(Region = rownames(table2)) %>%
    select(Region, everything())
  
  list(
    has_prot2 = has_p2,
    raw = list(dna = dna, brdu = brdu, prot1 = p1, prot2 = if (has_p2) p2 else NULL),
    table1 = table1,
    table2 = table2
  )
}

# ---- Tab 1 data: lane-colored alpha, 0 => white ----
fiber_plot_df_alpha <- function(raw, has_prot2,
                               dna_thres, brdu_thres, p1_thres, p2_thres,
                               label_dna, label_brdu, label_p1, label_p2,
                               col_dna, col_brdu, col_p1, col_p2) {
  
  n <- length(raw$dna)
  pos <- seq_len(n)
  
  # mask below threshold to 0 (gap)
  dna  <- ifelse(raw$dna  >= dna_thres,  raw$dna,  0)
  brdu <- ifelse(raw$brdu >= brdu_thres, raw$brdu, 0)
  p1   <- ifelse(raw$prot1 >= p1_thres,  raw$prot1, 0)
  
  lane_names <- c(label_dna, label_brdu, label_p1)
  lane_colors <- c(col_dna, col_brdu, col_p1)
  
  lanes <- list(dna = dna, brdu = brdu, p1 = p1)
  
  if (has_prot2) {
    p2 <- ifelse(raw$prot2 >= p2_thres, raw$prot2, 0)
    lanes$p2 <- p2
    lane_names <- c(lane_names, label_p2)
    lane_colors <- c(lane_colors, col_p2)
  }
  
  df <- bind_rows(lapply(seq_along(lanes), function(i) {
    vals <- lanes[[i]]
    nm <- lane_names[i]
    tibble(position = pos, lane = nm, value = vals, lane_color = lane_colors[i])
  })) %>%
    group_by(lane) %>%
    mutate(alpha = ifelse(value <= 0, 0, value / max(value, na.rm = TRUE))) %>%
    ungroup() %>%
    mutate(lane = factor(lane, levels = lane_names))
  
  df
}

plot_fiber_lanes_alpha <- function(df) {
  ggplot(df, aes(x = position, y = 1)) +
    geom_raster(aes(alpha = alpha), fill = df$lane_color) +
    facet_grid(lane ~ ., scales = "free_y", switch = "y") +
    scale_y_continuous(NULL, breaks = NULL) +
    scale_alpha(range = c(0, 1), guide = "none") +
    labs(x = "Pixel position along fiber") +
    theme_minimal(base_size = 12) +
    theme(
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, margin = margin(r = 0)), # tighter
      strip.background = element_blank(),
      panel.spacing.y = unit(0.10, "lines"),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 2),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA)
    )
}

plot_table2_bars <- function(table2) {
  metric_cols <- intersect(c("Prot1Percent", "Prot2Percent", "ColocalizationPercent"), names(table2))
  
  df <- table2 %>%
    filter(Region %in% c("Replicated", "Forks", "Unreplicated")) %>%
    select(Region, all_of(metric_cols)) %>%
    pivot_longer(-Region, names_to = "Metric", values_to = "Percent")
  
  ggplot(df, aes(x = Region, y = Percent, fill = Metric)) +
    geom_col(position = "dodge") +
    scale_fill_viridis_d(option = "D", end = 0.9) +
    labs(x = NULL, y = "Percent of protein signal", fill = NULL) +
    theme_minimal(base_size = 12)
}

ui <- fluidPage(
  theme = bslib::bs_theme(version = 5),
  uiOutput("app_title"),
  
  fluidRow(
    column(
      4,
      bslib::card(
        bslib::card_body(
          textInput("expt", "Experiment name", value = "MyExperiment"),
          fileInput("fiberfile", "Upload fiber .txt", accept = c(".txt", ".tsv", ".tab")),
          verbatimTextOutput("detected_channels", placeholder = TRUE),
          
          bslib::accordion(
            bslib::accordion_panel("Channel mapping", uiOutput("mapping_ui")),
            bslib::accordion_panel("Thresholds (analysis)",
              numericInput("brdu_thres", "BrdU threshold (tract calling)", value = 200, min = 0),
              numericInput("p1_thres", "Protein 1 threshold", value = 200, min = 0),
              numericInput("p2_thres", "Protein 2 threshold", value = 200, min = 0)
            ),
            bslib::accordion_panel("Thresholds (display)",
              numericInput("dna_thres_disp", "DNA display threshold", value = 0, min = 0)
            ),
            bslib::accordion_panel("Fork + smoothing",
              numericInput("pR", "Fork pixels into replicated region (pR)", value = 10, min = 0),
              numericInput("pU", "Fork pixels into unreplicated region (pU)", value = 10, min = 0),
              numericInput("smooth", "Smooth it (close gaps <= this many pixels)", value = 0, min = 0)
            ),
            bslib::accordion_panel("Tab 1 edits",
              textInput("col_dna",  "DNA color (hex)", value = "#1f77b4"),
              textInput("col_brdu", "BrdU color (hex)", value = "#2ca02c"),
              textInput("col_p1",   "Protein 1 color (hex)", value = "#d62728"),
              textInput("col_p2",   "Protein 2 color (hex)", value = "#ff7f0e"),
              uiOutput("color_swatches"),
              actionButton("rename_lanes", "Rename lanes"),
              helpText("These settings affect Tab 1 only.")
            ),
            bslib::accordion_panel("Export",
              downloadButton("download_pdf_tab1", "Save Tab 1 as PDF"),
              downloadButton("download_pdf_tab2", "Save Tab 2 as PDF")
            ),
            open = "Channel mapping"
          ),

          # Keep lane label values in hidden inputs; edited via modal dialog
          tags$div(style="display:none;",
            textInput("label_dna",  "DNA lane label", value = "DNA"),
            textInput("label_brdu", "BrdU lane label", value = "BrdU"),
            textInput("label_p1",   "Protein 1 lane label", value = "Protein 1"),
            textInput("label_p2",   "Protein 2 lane label", value = "Protein 2")
          )

        )
      )
    ),
    
    column(
      8,
      tabsetPanel(
        tabPanel(
          "Tab 1: Fiber lanes",
          br(),
          fluidRow(
            column(8),
            column(4, sliderInput("px_per_pos", "Horizontal zoom (px per fiber pixel)", min = 1, max = 8, value = 2, step = 1))
          ),
          tags$div(
            style = "overflow-x: auto; overflow-y: hidden; border: 1px solid #eee; padding: 6px;",
            uiOutput("fiber_plot_ui")
          ),
          br(),
          h4("table1 (first 25 rows)"),
          tableOutput("table1_head")
        ),
        tabPanel(
          "Tab 2: Region summary",
          br(),
          plotOutput("table2_plot", height = "450px"),
          br(),
          h4("table2"),
          tableOutput("table2_out")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  
  output$app_title <- renderUI({
    titlePanel(if (!is.null(input$expt) && nzchar(input$expt)) input$expt else "ODD-BLOBS Fiber Visualizer")
  })
  
  fiber_df <- reactiveVal(NULL)
  
  observeEvent(input$fiberfile, {
    req(input$fiberfile$datapath)
    df <- read_fiber_table(input$fiberfile$datapath)
    fiber_df(df)
  })
  
  output$mapping_ui <- renderUI({
    df <- fiber_df()
    if (is.null(df)) return(helpText("Upload a file to choose channels."))
    
    channel_cols <- which(grepl("^Channel", names(df)))
    if (length(channel_cols) < 3) return(helpText("No 'Channel ...' columns detected."))
    
    choices <- setNames(seq_along(channel_cols), paste0("Channel ", seq_along(channel_cols), " (", names(df)[channel_cols], ")"))
    
    tagList(
      selectInput("dna_idx",  "DNA lane uses:",  choices = choices, selected = 1),
      selectInput("brdu_idx", "BrdU lane uses:", choices = choices, selected = min(3, length(channel_cols))),
      selectInput("p1_idx",   "Protein 1 lane uses:", choices = choices, selected = min(2, length(channel_cols))),
      selectInput("p2_idx",   "Protein 2 lane uses (optional):", choices = c("None" = NA, choices), selected = if (length(channel_cols) >= 4) 4 else NA)
    )
  })
  
  output$detected_channels <- renderText({
    df <- fiber_df()
    if (is.null(df)) return("Upload a file to detect channels.")
    channel_cols <- which(grepl("^Channel", names(df)))
    paste0("Detected Channel columns: ", length(channel_cols), "\n",
           paste(names(df)[channel_cols], collapse = "\n"))
  })
  
  
  output$color_swatches <- renderUI({
    # Live preview of the hex inputs (Tab 1 only)
    swatch <- function(col) {
      tags$div(style = paste0("width:18px;height:18px;display:inline-block;border:1px solid #ccc;background:", col, ";"))
    }
    col_dna  <- normalize_hex(input$col_dna,  "#1f77b4")
    col_brdu <- normalize_hex(input$col_brdu, "#2ca02c")
    col_p1   <- normalize_hex(input$col_p1,   "#d62728")
    col_p2   <- normalize_hex(input$col_p2,   "#ff7f0e")

    tags$div(
      style="display:flex; gap:10px; align-items:center; flex-wrap:wrap; margin-top:6px;",
      tags$span("Preview:"),
      swatch(col_dna),  tags$span(style="font-size:12px;", "DNA"),
      swatch(col_brdu), tags$span(style="font-size:12px;", "BrdU"),
      swatch(col_p1),   tags$span(style="font-size:12px;", "Prot1"),
      swatch(col_p2),   tags$span(style="font-size:12px;", "Prot2")
    )
  })

  observeEvent(input$rename_lanes, {
    showModal(modalDialog(
      title = "Rename lanes (Tab 1)",
      textInput("label_dna_modal",  "DNA lane label",  value = input$label_dna  %||% "DNA"),
      textInput("label_brdu_modal", "BrdU lane label", value = input$label_brdu %||% "BrdU"),
      textInput("label_p1_modal",   "Protein 1 lane label", value = input$label_p1 %||% "Protein 1"),
      textInput("label_p2_modal",   "Protein 2 lane label", value = input$label_p2 %||% "Protein 2"),
      footer = tagList(
        actionButton("save_lane_labels", "Save", class = "btn-primary"),
        modalButton("Cancel")
      ),
      easyClose = TRUE
    ))
  })

  observeEvent(input$save_lane_labels, {
    updateTextInput(session, "label_dna",  value = input$label_dna_modal)
    updateTextInput(session, "label_brdu", value = input$label_brdu_modal)
    updateTextInput(session, "label_p1",   value = input$label_p1_modal)
    updateTextInput(session, "label_p2",   value = input$label_p2_modal)
    removeModal()
  })

# ---- AUTO-RUN: recompute whenever inputs change ----
  results <- reactive({
    df <- fiber_df()
    req(df)
    req(input$dna_idx, input$brdu_idx, input$p1_idx) # mapping must exist
    
    run_oddblobs_with_mapping(
      df = df,
      dna_idx = as.integer(input$dna_idx),
      brdu_idx = as.integer(input$brdu_idx),
      p1_idx = as.integer(input$p1_idx),
      p2_idx = if (is.na(input$p2_idx)) NA_integer_ else as.integer(input$p2_idx),
      brdu_thres = input$brdu_thres,
      p1_thres = input$p1_thres,
      p2_thres = input$p2_thres,
      pR = input$pR,
      pU = input$pU,
      smooth_val = input$smooth
    )
  })
  
  # dynamic scrollable plot width
  output$fiber_plot_ui <- renderUI({
    res <- results()
    req(res)
    n <- length(res$raw$dna)
    w <- max(900, n * input$px_per_pos)
    plotOutput("fiber_plot_static", height = "320px", width = paste0(w, "px"))
  })
  
  output$fiber_plot_static <- renderPlot({
    res <- results()
    req(res)
    
    dfp <- fiber_plot_df_alpha(
      raw = res$raw,
      has_prot2 = res$has_prot2,
      dna_thres = input$dna_thres_disp,
      brdu_thres = input$brdu_thres,
      p1_thres = input$p1_thres,
      p2_thres = input$p2_thres,
      label_dna = input$label_dna,
      label_brdu = input$label_brdu,
      label_p1 = input$label_p1,
      label_p2 = input$label_p2,
      col_dna = normalize_hex(input$col_dna, "#1f77b4"),
      col_brdu = normalize_hex(input$col_brdu, "#2ca02c"),
      col_p1 = normalize_hex(input$col_p1, "#d62728"),
      col_p2 = normalize_hex(input$col_p2, "#ff7f0e")
      )
    
    plot_fiber_lanes_alpha(dfp) + ggtitle(input$expt)
  }, res = 96)
  
  output$table1_head <- renderTable({
    res <- results()
    req(res)
    head(res$table1, 25)
  })
  
  output$table2_plot <- renderPlot({
    res <- results()
    req(res)
    plot_table2_bars(res$table2) + ggtitle(input$expt)
  })
  
  output$table2_out <- renderTable({
    res <- results()
    req(res)
    res$table2
  })
  
  # ---- Download PDF of Tab 1 ----
  output$download_pdf_tab1 <- downloadHandler(
    filename = function() {
      paste0(input$expt, "_fiber_plot.pdf")
    },
    content = function(file) {
      res <- results()
      req(res)
      
      dfp <- fiber_plot_df_alpha(
        raw = res$raw,
        has_prot2 = res$has_prot2,
        dna_thres = input$dna_thres_disp,
        brdu_thres = input$brdu_thres,
        p1_thres = input$p1_thres,
        p2_thres = input$p2_thres,
        label_dna = input$label_dna,
        label_brdu = input$label_brdu,
        label_p1 = input$label_p1,
        label_p2 = input$label_p2,
      col_dna = normalize_hex(input$col_dna, "#1f77b4"),
      col_brdu = normalize_hex(input$col_brdu, "#2ca02c"),
      col_p1 = normalize_hex(input$col_p1, "#d62728"),
      col_p2 = normalize_hex(input$col_p2, "#ff7f0e")
      )
      
      p <- plot_fiber_lanes_alpha(dfp) + ggtitle(input$expt)
      
      # Use a wide PDF so it prints nicely; scale with fiber length a bit.
      n <- length(res$raw$dna)
      width_in <- max(11, min(60, n / 150))  # cap so it doesn't get absurd
      height_in <- if (res$has_prot2) 4 else 3
      
      pdf(file, width = width_in, height = height_in, useDingbats = FALSE)
      print(p)
      dev.off()
    }
  )
  
  output$download_pdf_tab2 <- downloadHandler(
    filename = function() paste0(input$expt, "_region_summary.pdf"),
    content = function(file) {
      res <- results()
      req(res)
      
      p <- plot_table2_bars(res$table2) + ggtitle(input$expt)
      
      pdf(file, width = 8.5, height = 5.0, useDingbats = FALSE)
      print(p)
      dev.off()
    }
  )
}

shinyApp(ui, server)