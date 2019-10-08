library(shiny)

shinyServer(function(session, input, output) {

    set_logging_session()

    openTab <- reactiveVal()
    observe({ openTab(input$navbar) })
    
    log_message("Setting up table module")
    table_vars <- callModule(module=table_panel, id="Table", datasets=datasets, open_tab=openTab)
    callModule(module=table_vis, id="Table", table_vars=table_vars)
    
    log_message("Setting up enrichment module")
    enrich_vars <- callModule(enrichment_server, id="Enrichment", table_vars=table_vars)
    callModule(
        module=enrichment_panel,
        id="Enrichment",
        dataset=input$dataset,
        enrich_vals=enrich_vars
    )

    log_message("Setting up MultivarVis module")
    callModule(module=pca_panel, id="MultivarVis", table_vars=table_vars)
    
    log_message("Setting up SingleFeature module")
    callModule(module=single_feature_panel, id="SingleFeature", table_vars=table_vars, datasets=datasets)

    log_message("Setting up navbar")
    observeEvent(openTab(), {
        updateTabsetPanel(session, "navbar", openTab())
    })
})
