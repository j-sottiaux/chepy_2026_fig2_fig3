# Comprehensive Omics Data Enrichment eXplorer (CODEX) — Shiny app

# Load required libraries ----
{
  library(tidyverse)
  library(ggplot2)
  library(plotly)
  library(ggrepel)

  library(shiny)
  library(shinyWidgets)
  library(bslib)
}

gc()
set.seed(12345)

# Source functions / objects
source(file = "scripts/00_functions.R")

app_version <- "v1.3.1"

# Theme ----
theme_codex <- bs_theme(
  version = 5,
  bootswatch = "litera",
  bg = "#FFFFFF",
  fg = "#111827",
  primary = "#2563EB",
  base_font = font_google("Inter")
) |>
  bs_add_rules("
    /* App chrome */
    .codex-header { padding: 12px 0 6px 0; }
    .codex-footer { border-top: 1px solid #E5E7EB; padding: 10px 0; opacity: .8; font-size: 12px; }
    .card { border-color: #E5E7EB; box-shadow: none; }

    /* Header harmonization */
    .card-header,
    .accordion-button{
      background-color: rgba(var(--bs-secondary-rgb), 0.08);
      border-bottom: 1px solid #E5E7EB;
    }
    .accordion-button:not(.collapsed){
      background-color: rgba(var(--bs-secondary-rgb), 0.12);
      color: inherit;
    }

    /* Sidebar info card look */
    .card.sidebar-info {
      border: 1px solid rgba(var(--bs-info-rgb), 0.55) !important;
    }
    .card.sidebar-info .card-header {
      background-color: rgba(var(--bs-info-rgb), 0.08) !important;
    }

    /* Sidebar header cohesion */
    .sidebar-head-title { font-weight: 700; font-size: 18px; letter-spacing: .02em; line-height: 1.2; }
    .sidebar-head-sub   { font-size: 12px; opacity: .75; margin-top: 2px; line-height: 1.2; }

    /* Narrower page margins */
    .container, .container-fluid { padding-left: .5rem; padding-right: .5rem; }

    /* Tighter spacing */
    .bslib-sidebar-layout { gap: .75rem; }

    /* Volcano action row (top) */
    .plot-actions-top{
      display:flex;
      justify-content:space-between;
      align-items:center;
      gap:.5rem;
      flex-wrap:wrap;
      margin-bottom:.5rem;
    }
    .plot-actions-left{ display:flex; gap:.5rem; align-items:center; flex-wrap:wrap; }
    .plot-actions-right{ display:flex; gap:.75rem; align-items:center; flex-wrap:wrap; }

    /* Enrichment action row (TOP-LEFT) */
    .enrich-actions-top{
      display:flex;
      justify-content:flex-start;
      align-items:center;
      gap:.5rem;
      flex-wrap:wrap;
      margin-bottom:.5rem;
    }

    /* GeneCards row */
    .genecards-row{ display:flex; gap:.5rem; align-items:center; flex-wrap:wrap; }

    /* Scrollable main area */
    #mainScroll{
      overflow-y: auto;
      height: calc(100vh - 160px);
      padding-bottom: .5rem;
    }
  ")


# UI ----
ui <- page_fillable(
  theme = theme_codex,
  tags$div(
    class = "codex-header",
    tags$div(
      style = "display:flex; align-items:center; gap:12px;",
      tags$div(
        tags$div("CODEX", style = "font-size:30px; font-weight:700; line-height:1.05;"),
        tags$div(
          "Comprehensive Omics Data Enrichment eXplorer",
          style = "font-size:15px; opacity:0.75; line-height:1.1;"
        )
      )
    )
  ),

  # JS handler to scroll main panel to top when volcano opens
  tags$script(HTML("
    Shiny.addCustomMessageHandler('scrollMainTop', function(x){
      var el = document.getElementById('mainScroll');
      if(el) el.scrollTo({top: 0, behavior: 'smooth'});
    });
  ")),
  layout_sidebar(
    fillable = TRUE,
    sidebar = sidebar(
      width = 360,
      card(
        class = "sidebar-info",
        card_header(
          tags$div(
            tags$div(class = "sidebar-head-title", "⚙️ PARAMETERS"),
            tags$div(class = "sidebar-head-sub", "4 criterias required to generate interactive plots.")
          )
        ),
        card_body(
          pickerInput(
            "db_name",
            tags$strong("1️⃣ Reference database"),
            choices = names(ref_proteome)
          ),
          hr(),
          pickerInput(
            "bio_cond",
            tags$strong("2️⃣ Biological condition"),
            choices = c("dcSSc_ATAp", "dcSSc_ATAn", "lcSSc_ACA")
          ),
          hr(),
          pickerInput(
            "cell_comp",
            tags$strong("3️⃣ Cellular compartment"),
            choices = c("nucleus", "cytoplasm")
          ),
          hr(),
          pickerInput(
            "enrich_cat",
            tags$strong("4️⃣ Core enrichment genes"),
            choices = c("gsea", "ora", "gsea & ora")
          ),
          div(style = "height:.5rem;"),
          actionButton("refresh_btn", "Generate plot", class = "btn btn-primary w-100")
        )
      )
    ),
    tags$div(
      id = "mainScroll",
      uiOutput("main_layout_ui")
    )
  ),
  tags$footer(
    class = "codex-footer",
    tags$div(
      style = "display:flex; justify-content:space-between; gap:12px; flex-wrap:wrap;",
      tags$span(paste("CODEX", app_version)),
      tags$span(paste0(
        "© ", format(Sys.Date(), "%Y"),
        " • J.Sottiaux • Institute for Translational Research in Inflammation (INFINITE) - U1286 • All rights reserved."
      ))
    )
  )
)

# Server ----
server <- function(input, output, session) {
  padj_threshold <- 0.05

  rv <- reactiveValues(
    enrich_df = NULL,
    selected_pathway = NULL,
    enrich_gg = NULL,
    volcano_gg = NULL,
    gene_clicked = NULL,
    volcano_open = FALSE,
    params = NULL
  )

  observeEvent(input$refresh_btn, {
    rv$selected_pathway <- NULL
    rv$gene_clicked <- NULL
    rv$volcano_gg <- NULL
    rv$volcano_open <- FALSE

    rv$params <- list(
      enrich_cat = input$enrich_cat,
      db_name    = input$db_name,
      cell_comp  = input$cell_comp,
      bio_cond   = input$bio_cond
    )

    enrich_name <- paste0(rv$params$db_name, "_", rv$params$cell_comp, "_", rv$params$bio_cond)
    df <- master_enrich_data[[enrich_name]]
    validate(need(!is.null(df), "No enrichment table found for this selection."))

    thr <- -log10(padj_threshold)

    rv$enrich_df <- df %>%
      mutate(
        logGSEA = -log10(gsea_padj + 1e-10),
        logORA = -log10(ora_padj + 1e-10),
        co_signif = if_else(
          is.finite(logGSEA) & is.finite(logORA) & logGSEA > thr & logORA > thr,
          TRUE, FALSE
        ),
        pathway_relation = factor(pathway_relation, levels = names(pathways_colors))
      )

    rv$enrich_gg <- create_codex_enrichment(
      enrich_category     = rv$params$enrich_cat,
      db_name             = rv$params$db_name,
      cell_compartment    = rv$params$cell_comp,
      condition           = rv$params$bio_cond,
      pathway_description = "",
      padj_threshold      = padj_threshold,
      save_path           = NULL # no autosave
    ) +
      labs(title = NULL, subtitle = NULL)
  })

  output$main_layout_ui <- renderUI({
    # Landing state (before first run)
    if (is.null(rv$enrich_df) || is.null(rv$params) || is.null(rv$enrich_gg)) {
      return(
        card(
          card_header(tags$strong("CODEX")),
          card_body(
            tags$div(
              class = "alert alert-secondary",
              "Click “Generate plot” to load enrichment results and enable interactive exploration."
            )
          )
        )
      )
    }

    # Explicit heights so the plot always reserves space reliably
    enrich_h <- if (is.null(rv$selected_pathway)) "72vh" else "58vh"
    volcano_h <- "58vh"

    enrich_header <- tagList(
      tags$strong("Enrichment integration"),
      tags$span(
        style = "opacity:.75; font-weight:400; margin-left:8px;",
        paste0(
          toupper(rv$params$db_name), " / ", rv$params$cell_comp, " / ", rv$params$bio_cond,
          " / padj≤", padj_threshold
        )
      ),
      tags$span(
        style = "opacity:.65; font-weight:400; margin-left:8px;",
        "• click a co-significant dot to open volcano"
      )
    )

    enrich_card <- card(
      card_header(enrich_header),
      card_body(
        div(
          class = "enrich-actions-top",
          downloadButton("dl_enrich_png", "Save PNG", class = "btn btn-info btn-sm")
        ),
        plotlyOutput("enrich_plot", height = enrich_h)
      )
    )

    if (is.null(rv$selected_pathway)) {
      return(enrich_card)
    }

    volcano_title <- tagList(
      tags$strong("Volcano plot"),
      tags$span(style = "opacity:.75; font-weight:400; margin-left:8px;", rv$selected_pathway)
    )

    volcano_body <- tagList(
      div(
        class = "plot-actions-top",
        div(
          class = "plot-actions-left genecards-row",
          downloadButton("dl_volcano_png", "Save PNG", class = "btn btn-info btn-sm"),
          actionButton("clear_pathway", "Clear pathway", class = "btn btn-warning btn-sm"),
          actionButton("open_genecards_btn", "Open GeneCards", class = "btn btn-primary btn-sm"),
          uiOutput("selected_gene_ui")
        )
      ),
      plotlyOutput("volcano_plot", height = volcano_h)
    )

    tagList(
      accordion(
        id = "volcano_acc",
        open = if (isTRUE(rv$volcano_open)) "volcano_item" else character(0),
        accordion_panel(
          title = volcano_title,
          value = "volcano_item",
          volcano_body
        )
      ),
      enrich_card
    )
  })

  output$enrich_plot <- renderPlotly({
    req(rv$enrich_df)

    df <- rv$enrich_df
    thr <- -log10(padj_threshold)

    x_lim <- round(max(df$logGSEA, na.rm = TRUE) + 3)
    y_lim <- round(max(df$logORA, na.rm = TRUE) + 3)

    df_bg <- df %>%
      dplyr::filter(!co_signif) %>%
      dplyr::filter(is.finite(logGSEA), is.finite(logORA))

    df_sig <- df %>%
      dplyr::filter(co_signif) %>%
      dplyr::filter(is.finite(logGSEA), is.finite(logORA)) %>%
      dplyr::mutate(
        hover = paste0(
          "<b>", Description, "</b>",
          "<br><b>-log10 GSEA padj:</b> ", round(logGSEA, 3),
          "<br><b>-log10 ORA padj:</b> ", round(logORA, 3)
        )
      )

    pal <- pathways_colors[levels(df_sig$pathway_relation)]

    plotly::plot_ly(
      data = df_bg,
      x = ~logGSEA, y = ~logORA,
      type = "scatter", mode = "markers",
      source = "enrich",
      showlegend = FALSE,
      hoverinfo = "skip",
      marker = list(size = 6, color = "rgba(120,120,120,0.4)")
    ) %>%
      plotly::add_trace(
        data = df_sig,
        x = ~logGSEA, y = ~logORA,
        type = "scatter", mode = "markers",
        inherit = FALSE,
        color = ~pathway_relation,
        colors = pal,
        marker = list(size = 9, opacity = 0.9),
        text = ~hover, hoverinfo = "text",
        customdata = ~Description,
        showlegend = TRUE
      ) %>%
      plotly::layout(
        xaxis = list(title = "-log10(GSEA padj)", range = c(0, x_lim), zeroline = FALSE),
        yaxis = list(title = "-log10(ORA padj)", range = c(0, y_lim), zeroline = FALSE),
        legend = list(orientation = "h", x = 0, y = 1.08),
        shapes = list(list(
          type = "rect", layer = "below",
          x0 = thr, x1 = x_lim, y0 = thr, y1 = y_lim,
          line = list(color = "grey30", dash = "dash", width = 1),
          fillcolor = "rgba(200,200,200,0.25)"
        )),
        annotations = list(list(
          x = x_lim, y = thr - 0.2,
          text = "<b>co-significance area</b>",
          showarrow = FALSE, xanchor = "right",
          font = list(size = 10, color = "grey30")
        )),
        margin = list(l = 60, r = 20, t = 20, b = 55)
      ) %>%
      plotly::config(displayModeBar = TRUE)
  })

  observeEvent(plotly::event_data("plotly_click", source = "enrich"), {
    ev <- plotly::event_data("plotly_click", source = "enrich")
    req(ev, rv$params)

    clicked <- ev$customdata
    req(!is.null(clicked), nzchar(clicked))

    rv$selected_pathway <- clicked
    rv$gene_clicked <- NULL
    rv$volcano_open <- TRUE

    session$sendCustomMessage("scrollMainTop", list())

    rv$volcano_gg <- create_codex_volcano(
      enrich_category     = rv$params$enrich_cat,
      db_name             = rv$params$db_name,
      cell_compartment    = rv$params$cell_comp,
      condition           = rv$params$bio_cond,
      pathway_description = rv$selected_pathway,
      save_path           = NULL # no autosave
    ) +
      labs(title = NULL, subtitle = NULL)
  })

  observeEvent(input$clear_pathway, {
    rv$selected_pathway <- NULL
    rv$gene_clicked <- NULL
    rv$volcano_gg <- NULL
    rv$volcano_open <- FALSE
    session$sendCustomMessage("scrollMainTop", list())
  })

  output$volcano_plot <- renderPlotly({
    req(rv$volcano_gg)

    validate(need(exists("logFC_volcano_cutoff", inherits = TRUE), "logFC_volcano_cutoff not found."))
    validate(need(exists("padj_volcano_cutoff", inherits = TRUE), "padj_volcano_cutoff not found."))

    v1 <- -logFC_volcano_cutoff
    v2 <- logFC_volcano_cutoff
    h1 <- -log10(padj_volcano_cutoff)

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

    plotly::ggplotly(rv$volcano_gg, tooltip = "text", source = "volcano") %>%
      plotly::layout(shapes = shapes) %>%
      plotly::config(displayModeBar = TRUE)
  })

  observeEvent(plotly::event_data("plotly_click", source = "volcano"), {
    click_data <- plotly::event_data("plotly_click", source = "volcano")
    req(click_data)

    gene_id <- click_data$customdata
    req(!is.null(gene_id), nzchar(gene_id))

    rv$gene_clicked <- gene_id
  })

  output$selected_gene_ui <- renderUI({
    if (is.null(rv$gene_clicked)) {
      tags$span(style = "opacity:.7; font-size:12px;", "Click a volcano point to select a gene.")
    } else {
      tags$span(
        style = "font-size:12px;",
        tags$b("Selected: "),
        tags$code(rv$gene_clicked)
      )
    }
  })

  observeEvent(input$open_genecards_btn, {
    req(rv$gene_clicked)
    utils::browseURL(paste0(
      "https://www.genecards.org/cgi-bin/carddisp.pl?gene=",
      rv$gene_clicked
    ))
  })

  output$dl_enrich_png <- downloadHandler(
    filename = function() {
      req(rv$params)
      paste0(
        toupper(rv$params$db_name), "_enrichment_integration_",
        rv$params$bio_cond, "_", rv$params$cell_comp, ".png"
      )
    },
    content = function(file) {
      req(rv$enrich_gg)
      ggsave(file, plot = rv$enrich_gg, width = 8, height = 7.8, units = "in", dpi = 600)
    }
  )

  output$dl_volcano_png <- downloadHandler(
    filename = function() {
      req(rv$params, rv$selected_pathway)
      paste0(
        toupper(rv$params$db_name), "_", rv$params$enrich_cat, "_volcano_",
        rv$params$bio_cond, "_", str_replace_all(tolower(rv$selected_pathway), "\\s+", "_"),
        ".png"
      )
    },
    content = function(file) {
      req(rv$volcano_gg)
      ggsave(file, plot = rv$volcano_gg, width = 8, height = 7.8, units = "in", dpi = 600)
    }
  )
}

# Launch app ----
shinyApp(ui, server)
