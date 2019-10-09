table_stats_plots <- c(
    "PValueHistogram"
)

#' @import ggplot2
do_pval_hist <- function(data, stat_patterns, bins) {
    
    arg_col <- sprintf("%s.P.Value", stat_patterns[1])
    bel_col <- sprintf("%s.P.Value", stat_patterns[2])
    
    arg_phist <- ggplot(data, aes_string(x=arg_col)) + 
        geom_histogram(bins=bins) + 
        theme_classic() +
        ggtitle(sprintf("Condition: %s (filtered)", stat_patterns[1]))
    bel_phist <- ggplot(data, aes_string(x=bel_col)) + 
        geom_histogram(bins=bins) + 
        theme_classic() +
        ggtitle(sprintf("Condition: %s (filtered)", stat_patterns[2]))
    
    cowplot::plot_grid(arg_phist, bel_phist, ncol=2)
}


table_stats_panel_ui <- function(id) {
    ns <- NS(id)
    
    tabPanel(
        id,
        fluidPage(
            tags$head(
                tags$style(type="text/css", "select { max-width: 240px; }"),
                tags$style(type="text/css", ".span4 { max-width: 290px; }"),
                tags$style(type="text/css", ".well { max-width: 280px; }")
            ),
            div(
                style = "display:flex; align-items:flex-start",
                wellPanel(
                    style = "float:left;",
                    selectInput(ns("table_stats_plot"), "Plot type", choices=table_stats_plots, selected=table_stats_plots[1]),
                    numericInput(ns("table_stats_bins"), "Bins", value=50, min=1, max=200, step=5)
                ),
                fluidPage(
                    style = "flex-grow:1; resize:horizontal; overflow-x: hidden; overflow-y: hidden;",
                    fluidRow(uiOutput(ns("TableStatsOutput")))
                )
            )
        )
    ) 
} 

table_stats_panel <- function(input, output, session, table_vars) {
    
    output$TableStatsOutput <- renderUI({ 
        message("Table stats output ", input$table_stats_plot)
        plotOutput(session$ns(input$table_stats_plot)) 
    })
    output$PValueHistogram <- renderPlot({ 
        message("Trigger rendering")
        do_pval_hist(table_vars$cached_filtered_table(), table_vars$stat_base(), input$table_stats_bins)
    })
}


