library(shiny)
library(tidyverse)
library(PCAtools)

# source("CraftOmics/MultivarVis.R")

do_pca <- function(sdf, cond) {
    mv$pca(sdf, cond)
}

do_pca_scree <- function(sdf) {
    mv$plot_component_fraction(sdf, max_comps=9)
}

multi_vis_plots <- c("PCA", "Cluster", "Histogram")

multivarvis_panel_ui <- function(id, datasets, conditions) {
    ns <- NS(id)
    
    tabPanel(
        id,
        fluidPage(
            tags$head(
                tags$style(type="text/css", "select { max-width: 290px; }"),
                tags$style(type="text/css", ".span4 { max-width: 340px; }"),
                tags$style(type="text/css", ".well { max-width: 330px; }")
            ),
            div(
                style = "display:flex; align-items:flex-start",
                wellPanel(
                    style = "float:left;",
                    selectInput(ns("vis_type"), "Visualization", choices=multi_vis_plots, selected="PCA"),
                    selectInput(ns("pca_condition"), "Condition", choices=conditions, selected=conditions[1]),
                    selectInput(ns("custom_names"), "Custom names", choices=c("none", conditions), selected="none"),
                    selectInput(ns("legend_position"), "Legend position", choices=c("none", "right"), selected="none"),
                    selectInput(ns("omit_samples"), "Omit samples", choices=colnames(datasets[[1]]), selected = NULL, multiple = TRUE),
                    selectInput(ns("filter_type"), "Filter on type", choices=c("none", colnames(SummarizedExperiment::colData(datasets[[1]]))), selected="none"),
                    selectInput(ns("filter_type_levels"), "Levels to inspect", choices=NULL, selected=NULL, multiple=TRUE),
                    checkboxInput(ns("show_pca_settings"), "Show PCA settings", value=FALSE),
                    conditionalPanel(
                        sprintf("input['%s'] == 1", ns("show_pca_settings")),
                        selectInput(ns("shape_cond"), "Shape cond", choices=c("none", conditions), selected="none"),
                        numericInput(ns("var_filter"), "Variance filter", value=0.1, min=0, max=1, step=0.05),
                        selectInput(ns("comp1"), "Component 1", choices=paste0("PC", 1:9), selected="PC1"),
                        selectInput(ns("comp2"), "Component 2", choices=paste0("PC", 1:9), selected="PC2"),
                        checkboxInput(ns("pair_plot"), "Pairplot display", value=FALSE),
                        numericInput(ns("pair_plot_count"), "Pairplot count", value=5, min=1, max=20, step=1),
                        checkboxInput(ns("show_scree"), "Show Scree", value=FALSE),
                        numericInput(ns("pca_label_size"), "Label size", value=3),
                        numericInput(ns("pca_point_size"), "Point size", value=5)
                    ),
                    conditionalPanel(
                        sprintf("input['%s'] == 'Histogram'", ns("vis_type")),
                        numericInput(ns("table_stats_bins"), "Bins", value=50, min=1, max=200, step=5)
                    )
                ),
                fluidPage(
                    style = "flex-grow:1; resize:horizontal; overflow-x: hidden; overflow-y: hidden;",
                    conditionalPanel(
                        sprintf("input['%s'] == 'PCA'", ns("vis_type")),
                        conditionalPanel(
                            sprintf("input['%s'] == 0", ns("pair_plot")),
                            plotOutput(ns("PCAOutput"), height=800),
                            conditionalPanel(
                                sprintf("input['%s'] == 1", ns("show_scree")),
                                plotOutput(ns("PCAOutputScree"))
                            )
                        ),
                        conditionalPanel(
                            sprintf("input['%s'] == 1", ns("pair_plot")),
                            plotOutput(ns("PairPlot"), height=800)
                        )
                    ),
                    conditionalPanel(
                        sprintf("input['%s'] == 'Cluster'", ns("vis_type")),
                        plotOutput(ns("Dendogram"), height=800)
                    ),
                    conditionalPanel(
                        sprintf("input['%s'] == 'Histogram'", ns("vis_type")),
                        plotOutput(ns("Histogram"))
                    )
                )
            )
        )
    )
}

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

#' @import ggdendro ggplot2
#' @importFrom stats complete.cases
do_dendogram = function(raw_data_m, raw_color_levels, labels=NULL, pick_top_variance=NULL, title="Dendogram", omit_samples=NULL) {
    
    if (!is.null(omit_samples)) {
        data_m <- raw_data_m[, (!colnames(raw_data_m) %in% omit_samples)]
        color_levels <- raw_color_levels[colnames(raw_data_m) %in% colnames(data_m)]
    }
    else {
        data_m <- raw_data_m
        color_levels <- raw_color_levels
    }
    
    samples <- colnames(data_m)
    
    if (is.null(labels)) {
        labels <- samples
    }
    
    # Setup data
    expr_m_nona <- data_m[complete.cases(data_m),]
    
    # Calculate tree
    scaledTransposedMatrix <- scale(t(expr_m_nona), center=TRUE, scale=TRUE)
    hc <- stats::hclust(stats::dist(scaledTransposedMatrix), "ave")
    dhc <- stats::as.dendrogram(hc)
    # Note - Label order is shuffled within this object! Be careful with coloring.
    ddata <- ggdendro::dendro_data(dhc, type="rectangle")
    
    # Prepare for plotting
    cluster_label_order <- match(ddata$labels$label, samples)
    ddata$labels$color <- color_levels[cluster_label_order]
    ddata$labels$label <- labels[cluster_label_order]
    
    # Visualize
    plt <- ggplot(segment(ddata)) +
        geom_segment(aes(x=.data$x, y=.data$y, xend=.data$xend, yend=.data$yend)) +
        theme_dendro() +
        geom_text(data=label(ddata),
                  aes(x=.data$x, y=.data$y, label=.data$label, color=.data$color),
                  vjust=0.5, hjust=0, size=6) +
        coord_flip() +
        scale_y_reverse(expand=c(0.2, 0)) +
        scale_x_continuous(expand=c(0,1)) +
        ggtitle(title)
    plt
}

do_pca <- function(sdf, colinfo, color_col, var_filter=0.1, custom_names=NULL, shape_cond=NULL, 
                   legend_position="none", pc1="PC1", pc2="PC2", label_size=3, point_size=5) {
    
    if (is.null(custom_names)) {
        rownames(colinfo) <- colnames(sdf)
    }
    else {
        unique_names <- make.unique(custom_names)
        rownames(colinfo) <- unique_names
        colnames(sdf) <- unique_names
    }
    
    sdf <- sdf[complete.cases(sdf), ]
    p <- PCAtools::pca(sdf, metadata=colinfo, removeVar=var_filter)
    PCAtools::biplot(p, x=pc1, y=pc2, colby=color_col, shape=shape_cond, legendPosition = legend_position, 
                     pointSize = point_size, labSize = label_size)
}

do_pca_scree <- function(sdf, omit_samples=NULL, remove_var=0.1, components=1:10) {
    
    sdf <- sdf[complete.cases(sdf), ]
    p <- PCAtools::pca(sdf, removeVar=remove_var)
    PCAtools::screeplot(p, components=components)
}

do_pair_plot <- function(sdf, colinfo, color_col, pairs=5, omit_samples=NULL, var_filter=0.1) {
    
    rownames(colinfo) <- colnames(sdf)
    sdf <- sdf[complete.cases(sdf), ]
    p <- PCAtools::pca(sdf, metadata=colinfo, removeVar=var_filter)
    PCAtools::pairsplot(p, components = 1:pairs, colby=color_col, pointSize = 2)
}

multivarvis_panel <- function(input, output, session, table_vars, datasets) {
    
    observeEvent(input$mark_all_omitted, {
        updateSelectInput(session, "omit_samples", selected=colnames(datasets[[1]]))
    })
    
    observeEvent(input$filter_type, {
        if (input$filter_type == "none") {
            updateSelectInput(session, "filter_type_levels", selected=NULL, choices=NULL)
        }
        else {
            levels <- SummarizedExperiment::colData(datasets[[1]])[[input$filter_type]]
            updateSelectInput(session, "filter_type_levels", selected=levels, choices=levels)
        }
    })
    
    target_sdf <- reactive({
        sdf_raw <- table_vars$cached_sdf()
        
        if (input$filter_type != "none") {
            sdf_raw <- sdf_raw[, SummarizedExperiment::colData(table_vars$dataset())[[input$filter_type]] %in% input$filter_type_levels]
        }
        
        if (!is.null(input$omit_samples)) {
            sdf_raw[, (!colnames(sdf_raw) %in% input$omit_samples)]
        }
        else {
            sdf_raw
        }
    })
    
    target_ddf <- reactive({
        ddf_raw <- SummarizedExperiment::colData(table_vars$dataset()) %>% data.frame()
        if (any(colnames(table_vars$cached_sdf()) != colnames(target_sdf()))) {
            ddf_raw[colnames(table_vars$cached_sdf()) %in% colnames(target_sdf()), ]
        }
        else {
            ddf_raw
        }
    })
    
    output$Dendogram <- renderPlot({
        
        if (input$custom_names != "none") sample_names <- target_ddf()[[input$custom_names]]
        else sample_names <- NULL
        
        do_dendogram(
            target_sdf(),
            target_ddf() %>% dplyr::select(input$pca_condition) %>% unlist(),
            omit_samples=input$omit_samples,
            labels=sample_names
        )
    })
    
    output$PCAOutput <- renderPlot({
        
        if (input$custom_names != "none") sample_names <- target_ddf()[[input$custom_names]]
        else sample_names <- NULL

        if (input$shape_cond != "none") {
            shape <- input$shape_cond
        }
        else {
            shape <- NULL
        }
        
        do_pca(
            target_sdf(), 
            target_ddf(), 
            input$pca_condition, 
            custom_names=sample_names,
            shape_cond=shape,
            legend_position=input$legend_position,
            pc1=input$comp1, 
            pc2=input$comp2,
            var_filter=input$var_filter,
            label_size=input$pca_label_size,
            point_size=input$pca_point_size
        )
    })
    
    output$Histogram <- renderPlot({ 
        message("Trigger rendering")
        do_pval_hist(table_vars$cached_filtered_table(), table_vars$stat_base(), input$table_stats_bins)
    })
    
    output$PCAOutputScree <- renderPlot({
        do_pca_scree(target_sdf(), omit_samples=input$omit_samples)
    })
    
    output$PairPlot <- renderPlot({
        do_pair_plot(
            target_sdf(), 
            target_ddf(), 
            input$pca_condition,
            pairs=input$pair_plot_count,
            omit_samples=input$omit_samples,
            var_filter=input$var_filter
        )
    })
}





