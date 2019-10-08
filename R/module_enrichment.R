library(shiny)

goe_plot_types <- c(
    "GOdot",
    "GObar",
    "GOnet",
    "GOnetcirc",
    "GOheat",
    "GOmap",
    "GOupset",
    "GOplot"
)

gsea_plot_types <- c(
    "GOdot",
    "GOnet",
    "GOnetcirc",
    "GOheat",
    "GOmap",
    "GOridge",
    "GOgsea"
)

enrichment_panel_ui <- function(id) {
    
    ns <- NS(id)
    tabPanel(
        id,
        tags$head(tags$style(HTML(".shiny-split-layout > div { overflow: visible; }"))),
        tags$head(
            tags$style(type="text/css", "select { max-width: 240px; }"),
            tags$style(type="text/css", ".span4 { max-width: 290px; }"),
            tags$style(type="text/css", ".well { max-width: 280px; }")
        ),
        fluidPage(
            div(
                style = "display:flex; align-items:flex-start",
                wellPanel(
                    style = "float:left;",
                    actionButton(ns("perform_enrich"), "Perform enrichment"),
                    selectInput(ns("enrich_type"), "Enrichment type", choices = c("GOE", "GSEA"), selected = "GOE"),
                    selectInput(ns("ontology"), "Ontology", choices=c("MF", "CC", "BP", "ALL"), selected="MF"),
                    selectInput(ns("enrichment_plot"), "Enrichment plot", choices=goe_plot_types, selected=goe_plot_types[1]),
                    checkboxInput(ns("separate_fdr_display"), "Separate FDR display", value=FALSE),
                    conditionalPanel(
                        sprintf("input['%s'] == 1", ns("separate_fdr_display")),
                        sliderInput(ns("sep_fdr_cutoff"), "FDR cutoff", value=0.1, min=0, max=1, step=0.01)
                    ),
                    checkboxInput(ns("enrich_stat_options"), "Show stat options", value=FALSE),
                    conditionalPanel(
                        sprintf("input['%s'] == 1", ns("enrich_stat_options")),
                        checkboxInput(ns("enrich_fdr"), "Use BH FDR", value = TRUE),
                        sliderInput(ns("enrich_cutoff"), "Cutoff", min=0, max=1, step=0.01, value=0.1)
                    ),
                    checkboxInput(ns("enrich_display_options"), "Show display options", value=FALSE),
                    conditionalPanel(
                        sprintf("input['%s'] == 1", ns("enrich_display_options")),
                        selectInput(ns("displayed_plot"), "Displayed plot", choices=c("Both", "First", "Second"), selected="Both"),
                        numericInput(ns("plot_cols"), "Plot cols", value=1, step=1, min=1, max=2),
                        numericInput(ns("enrichment_plot_height"), "Plot height", value=500, min=0, max=5000, step=50),
                        numericInput(ns("max_display_terms"), "Display terms", value=30)
                    ),
                    div(textOutput(ns("enrichment_status")))
                ),
                fluidPage(
                    style = "flex-grow:1; resize:horizontal; overflow:hidden;",
                    fluidRow(
                        uiOutput(ns("EnrichmentPlot"))
                    )
                )
            )
        )
    )
} 

do_enrich_plot <- function(plot_func, first_enrich, second_enrich, plot_type, stat_bases, 
                           cols=2, fold_first=NULL, fold_second=NULL, ...) {
    
    titles <- stat_bases
    
    if (!(plot_type %in% c("Both", "First", "Second"))) {
        stop("Unknown plot type: ", plot_type)
    }
    
    make_plot <- function(enrich_obj, cutoff, title, plot_func, fold=NULL) {
        if (length(which(enrich_obj@result$p.adjust < cutoff)) == 0) {
            plt <- ggplot() + ggtitle(paste(title, "no enriched terms found!"))
        }
        else {
            if (!is.null(fold_first)) {
                plt <- plot_func(enrich_obj, foldChange=fold, ...)
            }
            else {
                plt <- plot_func(enrich_obj, ...)
            }
            plt <- plt + ggtitle(title)
        }
        plt
    }
    
    if (plot_type == "Both" || plot_type == "First") {
        if (!is.null(fold_first)) 
            plt1 <- make_plot(first_enrich, 0.1, titles[1], plot_func, fold=fold_first)
        else 
            plt1 <- make_plot(first_enrich, 0.1, titles[1], plot_func)
    }
    
    if (plot_type == "Both" || plot_type == "Second") {
        if (!is.null(fold_second)) 
            plt2 <- make_plot(second_enrich, 0.1, titles[2], plot_func, fold=fold_second)
        else 
            plt2 <- make_plot(second_enrich, 0.1, titles[2], plot_func)
    }
    
    # if (plot_type == "Both" || plot_type == "Second") {
    #     if (!is.null(fold_second)) {
    #         plt2 <- plot_func(second_enrich, foldChange=fold_second, ...)
    #     }
    #     else {
    #         plt2 <- plot_func(second_enrich, ...)
    #     }
    #     plt2 <- plt2 + ggtitle(titles[2])
    # }

    if (plot_type == "Both") {
        plot_grid(plt1, plt2, ncol=cols)
    }
    else if (plot_type == "First") {
        plt1
    }
    else if (plot_type == "Seconc") {
        plt2
    }
    else {
        stop("Unknown plot-type: ", plot_type)
    }
}

enrichment_server <- function(input, output, session, table_vars) {
    
    enrich_vals <- list()
    enrich_vals$out <- reactive({
        if (input$perform_enrich == 0) {
            return()
        }
        isolate({
            vars <- trigger_enrich(
                table_vars,
                table_vars$dataset(),
                input$enrich_type,
                input$enrich_fdr,
                input$enrich_cutoff,
                ontology=input$ontology,
                stat_patterns=table_vars$stat_base(),
                separate_fdr_display=input$separate_fdr_display,
                separate_fdr_cutoff=input$sep_fdr_cutoff
            )
            return(vars)
        })
    })
    
    return(enrich_vals)
}

trigger_enrich = function(table_vars, dataset, enrich_type, enrich_fdr, enrich_cutoff, ontology, stat_patterns, separate_fdr_display=FALSE, separate_fdr_cutoff=0.1) {

    enrich_vars <- list()
    enrich_vars$enrichment_status <- "Processing enrichment"
    
    if (!separate_fdr_display) {
        filter_df_first <- table_vars$cached_filtered_table()
        filter_df_second <- filter_df_first
    }
    else {
        filter_df_first <- table_vars$cached_full_table() %>% filter(UQ(as.name(sprintf("%s.adj.P.Val", stat_patterns[1]))) < separate_fdr_cutoff)
        filter_df_second <- table_vars$cached_full_table() %>% filter(UQ(as.name(sprintf("%s.adj.P.Val", stat_patterns[2]))) < separate_fdr_cutoff)
    }
    
    # Separate FDR display goes here
    
    showNotification(
        paste0("Performing ", enrich_type ," enrichment for ", nrow(filter_df_first), " and ", nrow(filter_df_second), 
               " features, it can take a few moments..."),
        duration = 4,
        type = "message"
    )
    
    enrich_vars$gene_list_first <- make_gene_list(filter_df_first, "ProteinID", sprintf("%s.logFC", stat_patterns[1]))
    enrich_vars$gene_list_second <- make_gene_list(filter_df_second, "ProteinID", sprintf("%s.logFC", stat_patterns[2]))
    universe <- make_universe(rowData(dataset) %>% data.frame(), "ProteinID")

    enrich_vars$enrichment_obj_first <- perform_enrichment(
        enrich_vars$gene_list_first, 
        universe, 
        enrich_type, 
        enrich_fdr, 
        enrich_cutoff,
        ontology=ontology
    )
    
    enrich_vars$enrichment_obj_second <- perform_enrichment(
        enrich_vars$gene_list_second, 
        universe, 
        enrich_type, 
        enrich_fdr, 
        enrich_cutoff,
        ontology=ontology
    )
    
    enrich_vars$enrichment_type <- enrich_type
    enrich_vars$enrichment_status <- paste(
        get_enrichment_status(enrich_type, enrich_vars$enrichment_obj_first, enrich_cutoff),
        get_enrichment_status(enrich_type, enrich_vars$enrichment_obj_second, enrich_cutoff),
        collapse="\n"
    )
    
    enrich_vars$stat_bases <- stat_patterns
    
    message("Enrichment done!")
    enrich_vars
}


enrichment_panel <- function(input, output, session, dataset, enrich_vals) {

    ns <- session$ns
    
    plotHeight <- reactive({
        input$enrichment_plot_height
    })

    output$EnrichmentPlot <- renderUI({
        
        if (is.null(enrich_vals$out())) {
            message("No enrichment assign, no rendering")
            message(ns("None"))
            plotOutput(ns("None"), height=plotHeight())
        }
        else {
            message("Going to render enrichment")
            message(ns(input$enrichment_plot))
            plotOutput(ns(input$enrichment_plot), height=plotHeight())
        }
    })

    output$None <- renderPlot({ ggplot() + ggtitle("Empty title") })
    
    output$GOdot <- renderPlot({
        
        do_enrich_plot(
            clusterProfiler::dotplot, 
            enrich_vals$out()$enrichment_obj_first, 
            enrich_vals$out()$enrichment_obj_second, 
            plot_type=input$displayed_plot, 
            cols=input$plot_cols,
            stat_bases=enrich_vals$out()$stat_bases,
            showCategory=input$max_display_terms
        )
    })
    
    output$GObar <- renderPlot({
        do_enrich_plot(
            barplot, 
            enrich_vals$out()$enrichment_obj_first, 
            enrich_vals$out()$enrichment_obj_second, 
            plot_type=input$displayed_plot, 
            cols=input$plot_cols,
            stat_bases=enrich_vals$out()$stat_bases,
            showCategory=input$max_display_terms
        )
    })
    
    output$GOnet <- renderPlot({
        
        do_enrich_plot(
            cnetplot,
            enrich_vals$out()$enrichment_obj_first, 
            enrich_vals$out()$enrichment_obj_second, 
            plot_type=input$displayed_plot,
            cols=input$plot_cols,
            fold_first=enrich_vals$out()$gene_list_first,
            fold_second=enrich_vals$out()$gene_list_second,
            stat_bases=enrich_vals$out()$stat_bases,
            showCategory=input$max_display_terms
        )
    })
    
    output$GOnetcirc <- renderPlot({
        do_enrich_plot(
            cnetplot,
            enrich_vals$out()$enrichment_obj_first, 
            enrich_vals$out()$enrichment_obj_second, 
            plot_type=input$displayed_plot,
            cols=input$plot_cols,
            fold_first=enrich_vals$out()$gene_list_first,
            fold_second=enrich_vals$out()$gene_list_second,
            circular=TRUE,
            colorEdge=TRUE,
            stat_bases=enrich_vals$out()$stat_bases,
            showCategory=input$max_display_terms
        )
    })
    
    output$GOheat <- renderPlot({
        do_enrich_plot(
            heatplot,
            enrich_vals$out()$enrichment_obj_first, 
            enrich_vals$out()$enrichment_obj_second, 
            plot_type=input$displayed_plot,
            cols=input$plot_cols,
            fold_first=enrich_vals$out()$gene_list_first,
            fold_second=enrich_vals$out()$gene_list_second,
            stat_bases=enrich_vals$out()$stat_bases,
            showCategory=input$max_display_terms
        )
    })
    
    output$GOmap <- renderPlot({
        do_enrich_plot(
            emapplot, 
            enrich_vals$out()$enrichment_obj_first, 
            enrich_vals$out()$enrichment_obj_second, 
            plot_type=input$displayed_plot, 
            cols=input$plot_cols,
            stat_bases=enrich_vals$out()$stat_bases,
            showCategory=input$max_display_terms
        )
    })
    
    output$GOridge <- renderPlot({
        do_enrich_plot(
            ridgeplot, 
            enrich_vals$out()$enrichment_obj_first, 
            enrich_vals$out()$enrichment_obj_second, 
            plot_type=input$displayed_plot, 
            cols=input$plot_cols,
            stat_bases=enrich_vals$out()$stat_bases,
            showCategory=input$max_display_terms
        )
    })
    
    output$GOgsea <- renderPlot({
        do_enrich_plot(
            gseaplot, 
            enrich_vals$out()$enrichment_obj_first, 
            enrich_vals$out()$enrichment_obj_second, 
            plot_type=input$displayed_plot, 
            cols=input$plot_cols, 
            geneSetID=1, 
            by="preranked",
            stat_bases=enrich_vals$out()$stat_bases,
            showCategory=input$max_display_terms
        )
    })
    
    output$GOplot <- renderPlot({
        
        do_enrich_plot(
            goplot, 
            enrich_vals$out()$enrichment_obj_first, 
            enrich_vals$out()$enrichment_obj_second, 
            plot_type=input$displayed_plot, 
            cols=input$plot_cols,
            stat_bases=enrich_vals$out()$stat_bases,
            showCategory=input$max_display_terms
        )
    })
    
    output$GOupset <- renderPlot({ upsetplot(enrich_vals$out()$enrichment_obj_first) + ggtitle(enrich_vals$out()$stat_bases[1]) })
    observeEvent(input$enrich_type, {

        if (is.null(input$enrich_type)) {
            message("NULL in enrichment")
        }
        else if (input$enrich_type == "GOE") {
            updateSelectInput(session, "enrichment_plot", choices=goe_plot_types)
        }
        else if (input$enrich_type == "GSEA") {
            updateSelectInput(session, "enrichment_plot", choices=gsea_plot_types)
        }
        else {
            message("Nothing to update to: ", vars$enrichment_type)
        }
    })
}


perform_enrichment <- function(gene_list, universe, enrich_type="GOE", do_fdr=TRUE, cutoff=0.1, ontology="MF") {
    
    if (do_fdr) {
        adjust_method="BH"
    }
    else {
        adjust_method="none"
    }
    
    key_type <- "TAIR"
    gene_set_permutations <- 1000
    bioc_db <- "org.At.tair.db"
    
    if (enrich_type == "GOE") {
        go_enrich <- clusterProfiler::enrichGO(
            names(gene_list), 
            OrgDb = bioc_db, 
            keyType = key_type, 
            universe = universe, 
            ont = ontology,
            pAdjustMethod = adjust_method, 
            pvalueCutoff = cutoff
        )
    }
    else if (enrich_type == "GSEA") {
        
        warning("Would need a way to reduce proteins before GSEA!")
        
        gene_list_wo_na <- gene_list[!is.na(gene_list)]
        if (length(gene_list_wo_na) != length(gene_list)) {
            message("ENRICHMENT: Removing ", length(gene_list) - length(gene_list_wo_na), " NA entries, reducing from ", 
                    length(gene_list), " to ", length(gene_list_wo_na))
        }
        message("ENRICHMENT: Number of non-unique: ", length(unique(names(gene_list))) - length(unique(names(gene_list_wo_na))))
        message("ENRICHMENT: Used adjustment method: ", adjust_method, " with cutoff: ", cutoff)
        
        go_enrich <- clusterProfiler::gseGO(
            gene_list_wo_na, 
            OrgDb = bioc_db, 
            nPerm = gene_set_permutations,
            keyType = key_type, 
            ont = ontology,
            pAdjustMethod = adjust_method, 
            pvalueCutoff = cutoff
        )
    }
    else {
        stop("Unsupported enrich_type: ", enrich_type)
    }
    
    go_enrich
}

get_enrichment_status <- function(enrich_type, enrich_obj, used_cutoff) {
    
    results_df <- enrich_obj@result
    sig_count <- results_df %>% filter(p.adjust < used_cutoff) %>% nrow()
    tot_count <- nrow(results_df)
    result_string <- paste0(enrich_type, " enriched terms: ", sig_count, " / ", tot_count, " Cutoff: ", used_cutoff)
    result_string
}

make_universe <- function(row_data, id_col) {
    all_ids <- unique(row_data[[id_col]])
    at_ids <- all_ids[grepl("^AT", all_ids)] %>% gsub("\\.\\d$", "", .)
    at_ids
}

make_gene_list <- function(row_data, id_col, stat_col, filter_col=NULL, filter_thres=0.1) {
    
    df <- row_data %>% data.frame()
    
    if (!is.null(filter_col)) {
        df <- df %>% filter(UQ(as.name(filter_col)) < filter_thres)
    }
    
    df <- df %>% dplyr::select(c(id_col, stat_col))
    colnames(df) <- c("id_col", "stat_col")
    df <- df %>%
        filter(grepl("^AT", id_col)) %>%
        arrange(desc(stat_col)) %>%
        mutate(id_col=gsub("\\.\\d$", "", id_col))
    geneList <- df$stat_col
    names(geneList) <- df$id_col
    geneList
}




