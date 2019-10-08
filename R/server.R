library(shiny)

server <- shinyServer(function(session, input, output) {

    global <- get_global()
    
    openTab <- reactiveVal()
    observe({ openTab(input$navbar) })
    
    table_vars <- callModule(module=table_panel, id="Table", datasets=global$datasets, open_tab=openTab)
    callModule(module=table_vis, id="Table", table_vars=table_vars)
    
    enrich_vars <- callModule(enrichment_server, id="Enrichment", table_vars=table_vars)
    callModule(
        module=enrichment_panel,
        id="Enrichment",
        dataset=input$dataset,
        enrich_vals=enrich_vars
    )

    callModule(module=pca_panel, id="MultivarVis", table_vars=table_vars)
    callModule(module=single_feature_panel, id="SingleFeature", table_vars=table_vars, datasets=global$datasets)

    observeEvent(openTab(), {
        updateTabsetPanel(session, "navbar", openTab())
    })
})
