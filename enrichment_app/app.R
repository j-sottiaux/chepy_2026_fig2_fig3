# Project setup  -----

# Load required libraries
library(shiny)
library(shinyWidgets)
library(bslib)
library(tidyverse)
library(ggplot2)
library(plotly)
library(ggrepel)

# Reproducibility check
gc()
set.seed(12345)

# Source functions
source(file = "scripts/01_functions.R")


# --- UI ----
ui <- fluidPage(
  
  titlePanel("Functional enrichment explorer"),
  
  sidebarLayout(
    sidebarPanel(
      h4("Select my filter options"),
      
      selectInput("enrich_cat", "Enrichment category", 
                  choices = c("gsea", "ora", "gsea & ora")),
      
      selectInput("db_name", "Reference database", 
                  choices = names(ref_proteome)),
      
      selectInput("cell_comp", "Cellular compartment", 
                  choices = c("nucleus", "cytoplasm")),
      
      selectInput("bio_cond", "Biological condition", 
                  choices = c("lcSSc_ACA", "dcSSc_ATAp", "dcSSc_ATAn")),
      
      selectInput("pathway_relation", "Pathway category", 
                  choices = c("DNA metabolism & processing", 
                              "Immune & stress reponse", 
                              "Other", 
                              "RNA metabolism & processing", 
                              "Translation & ribosome biology")),
      
      uiOutput("pathway_description_ui"),  # dynamically generated
      
      actionButton("refresh_btn", "Update my plot"),
      actionButton("save", "Save my plot")
    ),
    
    mainPanel(
      plotlyOutput("plot", height = "700px"),
      htmlOutput("gene_info")
    )
  )
)



# --- Server ----
server <- function(input, output, session) {
  
  # Dynamically populate "Pathway description" based on selected category
  output$pathway_description_ui <- renderUI({
    enrich_name <- paste0(input$db_name, "_", input$cell_comp, "_", input$bio_cond)
    df <- master_enrich_data[[enrich_name]]
    
    if (is.null(df)) return(NULL)
    
    filtered_df <- df %>% 
      filter(pathway_relation == input$pathway_relation)
    
    pathway_choices <- unique(filtered_df$Description)
    
    selectInput("pathway_desc", "Pathway description", 
                choices = pathway_choices,
                selected = pathway_choices[1])
  })
  
  # Dynamically update title
  output$plot_title <- renderText({
    paste("Volcano plot /", input$cell_comp, "_", input$bio_cond, "_HC")
  })
  
  # Reactive plot object
  plot_reactive <- reactive({
    req(input$pathway_desc)
    
    create_master_volcano(
      enrich_category = input$enrich_cat,
      db_name = input$db_name,
      cell_compartment = input$cell_comp,
      condition = input$bio_cond,
      pathway_description = input$pathway_desc
    )
  })
  
  # Render plotly plot
  output$plot <- renderPlotly({
    req(input$refresh_btn)
    isolate({
      ggplotly(plot_reactive(), tooltip = "text")
    })
  })
  
  # Save plot when button is clicked
  observeEvent(input$save, {
    create_master_volcano(
      enrich_category = input$enrich_cat,
      db_name = input$db_name,
      cell_compartment = input$cell_comp,
      condition = input$bio_cond,
      pathway_description = input$pathway_desc,
      save_path = "figures/"
    )
    showNotification("Plot saved", type = "message")
  })
  
  observeEvent(event_data("plotly_click"), {
    click_data <- event_data("plotly_click")
    if (is.null(click_data)) return()
    
    gene_id <- click_data$customdata  # no need to parse HTML
    
    gene_html <- paste0(
      "<b>Gene:</b> <a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=",
      gene_id, "' target='_blank'>", gene_id, "</a>"
    )
    
    output$gene_info <- renderUI({
      HTML(gene_html)
    })
  })
  
}

# --- Run app ----
shinyApp(ui = ui, server = server)