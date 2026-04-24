# Differential Analysis Explorer (DAE) — Shiny app

# Load required libraries ----
{
  library(tidyverse)
  library(ggplot2)
  library(plotly)

  library(shiny)
  library(shinyWidgets)
  library(bslib)

  library(ggrepel)
}

gc()
set.seed(12345)

# Source functions + data ----
source("scripts/01_functions.R")

conditions <- c("ATAp", "ATAn", "ACA", "HC", "NS")
condition_choices <- setNames(
  conditions,
  unname(proteo_cond_labels[conditions])
)

stopifnot(exists("proteo_toptables_extended", inherits = TRUE))
stopifnot(exists("create_dae_volcano", inherits = TRUE))
stopifnot(exists("logFC_threshold", inherits = TRUE))
stopifnot(exists("padj_threshold", inherits = TRUE))

app_version <- "v1.3.2"

# Theme ----
theme_dae <- bs_theme(
  version = 5,
  bootswatch = "litera",
  bg = "#FFFFFF",
  fg = "#111827",
  primary = "#2563EB",
  base_font = font_google("Inter")
) |>
  bs_add_rules("
    .dae-header { padding: 12px 0 6px 0; }
    .dae-footer { border-top: 1px solid #E5E7EB; padding: 10px 0; opacity: .8; font-size: 12px; }
    .card { border-color: #E5E7EB; box-shadow: none; }

    .card-header,
    .accordion-button{
      background-color: rgba(var(--bs-secondary-rgb), 0.08);
      border-bottom: 1px solid #E5E7EB;
    }
    .accordion-button:not(.collapsed){
      background-color: rgba(var(--bs-secondary-rgb), 0.12);
      color: inherit;
    }

    .card.sidebar-info {
      border: 1px solid rgba(var(--bs-info-rgb), 0.55) !important;
    }
    .card.sidebar-info .card-header {
      background-color: rgba(var(--bs-info-rgb), 0.08) !important;
    }

    .sidebar-head-title { font-weight: 700; font-size: 18px; letter-spacing: .02em; line-height: 1.2; }
    .sidebar-head-sub   { font-size: 12px; opacity: .75; margin-top: 2px; line-height: 1.2; }

    .container, .container-fluid { padding-left: .5rem; padding-right: .5rem; }
    .bslib-sidebar-layout { gap: .75rem; }

    .plot-actions-top{
      display:flex;
      justify-content:flex-start;
      align-items:center;
      gap:.5rem;
      flex-wrap:wrap;
      margin-bottom:.5rem;
    }

    .genecards-row { display:flex; gap:.5rem; align-items:center; flex-wrap:wrap; }

    /* KEY FIX: scroll container with a fixed viewport height */
    #mainScroll{
      overflow-y: auto;
      height: calc(100vh - 160px);
      padding-bottom: .5rem;
    }
  ")

# UI ----
ui <- page_fillable(
  theme = theme_dae,
  tags$div(
    class = "dae-header",
    tags$div(
      style = "display:flex; align-items:center; gap:12px;",
      tags$div(
        tags$div("Differential Analysis Explorer", style = "font-size:26px; font-weight:700; line-height:1.05;"),
        tags$div(
          "Explore pairwise contrasts and volcano plots from precomputed toptables with Limma (v3.66.0)",
          style = "font-size:14px; opacity:0.75; line-height:1.1;"
        )
      )
    )
  ),
  layout_sidebar(
    fillable = TRUE,
    sidebar = sidebar(
      width = 360,
      card(
        class = "sidebar-info",
        card_header(
          tags$div(
            tags$div(class = "sidebar-head-title", "⚙️ PARAMETERS"),
            tags$div(class = "sidebar-head-sub", "Choose compartment and contrast to display volcano plot.")
          )
        ),
        card_body(
          pickerInput(
            "cell_comp",
            tags$strong("1️⃣ Cellular compartment"),
            choices = c("nucleus", "cytoplasm"),
            selected = "nucleus"
          ),
          hr(),
          pickerInput(
            "cond_a",
            tags$strong("2️⃣ Condition A"),
            choices = condition_choices,
            selected = "ATAp"
          ),
          hr(),
          pickerInput(
            "cond_b",
            tags$strong("3️⃣ Condition B"),
            choices = condition_choices,
            selected = "HC"
          ),
          hr(),
          checkboxInput("use_shrunk", "Use logFC_shrunk", value = TRUE),
          div(style = "height:.5rem;"),
          actionButton("refresh_btn", "Generate volcano", class = "btn btn-primary w-100")
        )
      )
    ),
    tags$div(
      id = "mainScroll",
      uiOutput("main_layout_ui")
    )
  ),
  tags$footer(
    class = "dae-footer",
    tags$div(
      style = "display:flex; justify-content:space-between; gap:12px; flex-wrap:wrap;",
      tags$span(paste("Differential Analysis Explorer", app_version)),
      tags$span(paste0(
        "© ", format(Sys.Date(), "%Y"),
        "• J.Sottiaux • Institute for Translational Research in Inflammation (INFINITE) - U1286 • All rights reserved."
      ))
    )
  )
)

# Server ----
server <- function(input, output, session) {
  rv <- reactiveValues(
    volcano_plot = NULL,
    volcano_df = NULL,
    volcano_lfc_col = NULL,
    err = NULL,
    selected_gene = NULL
  )

  fetch_toptable_autoflip <- function(toptables, comp, a, b) {
    key_ab <- paste0(comp, "_", a, "_", b)
    key_ba <- paste0(comp, "_", b, "_", a)

    if (!is.null(toptables[[key_ab]])) {
      return(list(table = toptables[[key_ab]], key_ab = key_ab))
    }
    if (!is.null(toptables[[key_ba]])) {
      tt <- toptables[[key_ba]]
      if ("logFC" %in% names(tt)) tt$logFC <- -tt$logFC
      if ("logFC_shrunk" %in% names(tt)) tt$logFC_shrunk <- -tt$logFC_shrunk
      return(list(table = tt, key_ab = key_ab))
    }
    list(table = NULL, key_ab = key_ab)
  }

  observeEvent(input$refresh_btn, {
    rv$err <- NULL
    rv$selected_gene <- NULL

    tryCatch(
      {
        req(input$cell_comp, input$cond_a, input$cond_b)
        validate(need(input$cond_a != input$cond_b, "Condition A and Condition B must be different."))

        res <- fetch_toptable_autoflip(
          toptables = proteo_toptables_extended,
          comp = input$cell_comp,
          a = input$cond_a,
          b = input$cond_b
        )
        validate(need(!is.null(res$table), "No toptable found for this contrast in this compartment (direct or reverse)."))

        mini_list <- list()
        mini_list[[res$key_ab]] <- res$table

        out <- create_dae_volcano(
          toptables = mini_list,
          cell_compartment = input$cell_comp,
          condition_a = input$cond_a,
          condition_b = input$cond_b,
          save_path = NULL,
          use_shrunk = isTRUE(input$use_shrunk)
        )

        if (!is.list(out) || is.null(out$plot) || is.null(out$data)) {
          stop("create_dae_volcano() must return list(plot=<ggplot>, data=<data.frame>, lfc_col=<string>).")
        }

        needed_cols <- c("gene_id", "lfc", "mlog10", "diffexpressed", "hover_text", "gene_label")
        missing_cols <- setdiff(needed_cols, names(out$data))
        if (length(missing_cols) > 0) {
          stop("create_dae_volcano() returned data missing columns: ", paste(missing_cols, collapse = ", "))
        }

        out$data <- out$data %>% dplyr::mutate(.pt_id = dplyr::row_number())

        rv$volcano_plot <- out$plot
        rv$volcano_df <- out$data
        rv$volcano_lfc_col <- out$lfc_col
      },
      error = function(e) {
        rv$err <- conditionMessage(e)
        rv$volcano_plot <- NULL
        rv$volcano_df <- NULL
        rv$volcano_lfc_col <- NULL
        showNotification(paste("Volcano generation failed:", rv$err), type = "error", duration = NULL)
      }
    )
  })

  observeEvent(plotly::event_data("plotly_click", source = "volcano_src"), {
    ed <- plotly::event_data("plotly_click", source = "volcano_src")
    req(ed, rv$volcano_df)
    req(!is.null(ed$customdata), length(ed$customdata) >= 1)

    pt_id <- ed$customdata[[1]]
    req(length(pt_id) == 1, !is.na(pt_id))

    gene <- rv$volcano_df %>%
      dplyr::filter(.pt_id == pt_id) %>%
      dplyr::pull(gene_id)

    gene <- gene[1]
    if (!is.na(gene) && nzchar(gene)) rv$selected_gene <- gene
  })

  observeEvent(input$open_genecards, {
    req(rv$selected_gene)
    utils::browseURL(paste0(
      "https://www.genecards.org/cgi-bin/carddisp.pl?gene=",
      rv$selected_gene
    ))
  })

  output$selected_gene_ui <- renderUI({
    if (is.null(rv$selected_gene)) {
      tags$span(style = "opacity:.7; font-size:12px;", "Click a point to select a gene.")
    } else {
      tags$span(style = "font-size:12px;", tags$b("Selected: "), tags$code(rv$selected_gene))
    }
  })

  output$main_layout_ui <- renderUI({
    req(input$cell_comp, input$cond_a, input$cond_b)

    header <- tagList(
      tags$strong("Volcano plot"),
      tags$span(
        style = "opacity:.75; font-weight:400; margin-left:8px;",
        paste0(
          toupper(input$cell_comp), " / ",
          format_proteo_condition(input$cond_a),
          " vs ",
          format_proteo_condition(input$cond_b)
        )
      )
    )

    if (!is.null(rv$err)) {
      return(
        card(
          card_header(header),
          card_body(
            tags$div(
              class = "alert alert-danger",
              tags$b("Error while generating volcano: "),
              tags$code(rv$err)
            )
          )
        )
      )
    }

    if (is.null(rv$volcano_df) || is.null(rv$volcano_plot)) {
      return(
        card(
          card_header(header),
          card_body(
            tags$div(
              class = "alert alert-secondary",
              "Click “Generate volcano” to load a contrast and render the plot."
            )
          )
        )
      )
    }

    card(
      card_header(header),
      card_body(
        div(
          class = "plot-actions-top genecards-row",
          downloadButton("dl_volcano_png", "Save PNG", class = "btn btn-primary btn-sm"),
          actionButton("open_genecards", "Open GeneCards", class = "btn btn-primary btn-sm"),
          uiOutput("selected_gene_ui")
        ),
        # KEY FIX: explicit height (reserves space)
        plotlyOutput("volcano_plotly", height = "72vh")
      )
    )
  })

  output$volcano_plotly <- renderPlotly({
    req(rv$volcano_df)

    df <- rv$volcano_df
    lfc_col <- if (!is.null(rv$volcano_lfc_col)) rv$volcano_lfc_col else "logFC"

    df_bg <- df %>% dplyr::filter(diffexpressed == "no")
    df_up <- df %>% dplyr::filter(diffexpressed == "up")
    df_dn <- df %>% dplyr::filter(diffexpressed == "down")
    df_lab <- df %>% dplyr::filter(!is.na(gene_label))

    x_lim <- max(abs(df$lfc), na.rm = TRUE)
    x_lim <- if (is.finite(x_lim)) round(x_lim + 0.5, 1) else 2

    y_max <- max(df$mlog10, na.rm = TRUE)
    y_max <- if (is.finite(y_max)) ceiling(y_max) else 10

    v1 <- -logFC_threshold
    v2 <- logFC_threshold
    h1 <- -log10(padj_threshold)

    shapes <- list(
      list(
        type = "line", xref = "x", yref = "paper", x0 = v1, x1 = v1, y0 = 0, y1 = 1,
        line = list(color = "#dd9d6b", width = 1, dash = "dash")
      ),
      list(
        type = "line", xref = "x", yref = "paper", x0 = v2, x1 = v2, y0 = 0, y1 = 1,
        line = list(color = "#dd9d6b", width = 1, dash = "dash")
      ),
      list(
        type = "line", xref = "paper", yref = "y", x0 = 0, x1 = 1, y0 = h1, y1 = h1,
        line = list(color = "#dd9d6b", width = 1, dash = "dash")
      )
    )

    p <- plotly::plot_ly() %>%
      plotly::add_trace(
        data = df_bg,
        x = ~lfc, y = ~mlog10,
        type = "scatter", mode = "markers",
        showlegend = FALSE,
        hoverinfo = "skip",
        customdata = ~.pt_id,
        marker = list(size = 6, color = "rgba(170,170,170,0.4)")
      ) %>%
      plotly::add_trace(
        data = df_dn,
        x = ~lfc, y = ~mlog10,
        type = "scatter", mode = "markers",
        showlegend = FALSE,
        text = ~hover_text, hoverinfo = "text",
        customdata = ~.pt_id,
        marker = list(size = 8, color = "rgba(24,147,146,0.9)")
      ) %>%
      plotly::add_trace(
        data = df_up,
        x = ~lfc, y = ~mlog10,
        type = "scatter", mode = "markers",
        showlegend = FALSE,
        text = ~hover_text, hoverinfo = "text",
        customdata = ~.pt_id,
        marker = list(size = 8, color = "rgba(196,58,80,0.9)")
      ) %>%
      plotly::add_text(
        data = df_lab,
        x = ~lfc, y = ~mlog10,
        text = ~gene_label,
        textposition = "top center",
        showlegend = FALSE
      ) %>%
      plotly::layout(
        shapes = shapes,
        xaxis = list(title = lfc_col, range = c(-x_lim, x_lim), zeroline = FALSE),
        yaxis = list(title = "-log10(adj.P.Val)", range = c(0, y_max), zeroline = FALSE),
        margin = list(l = 60, r = 20, t = 20, b = 55)
      ) %>%
      plotly::config(displayModeBar = TRUE)

    p$x$source <- "volcano_src"
    p <- plotly::event_register(p, "plotly_click")
    p
  })

  output$dl_volcano_png <- downloadHandler(
    filename = function() {
      paste0(
        "DAE_",
        input$cell_comp, "_",
        safe_proteo_condition(input$cond_a), "_vs_",
        safe_proteo_condition(input$cond_b),
        ".png"
      )
    },
    content = function(file) {
      req(rv$volcano_plot)
      ggplot2::ggsave(file, plot = rv$volcano_plot, width = 8, height = 7.8, units = "in", dpi = 600)
    }
  )
}

shinyApp(ui, server)
