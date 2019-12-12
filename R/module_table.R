table_panel_ui <- function(
    id, timepoints, annot_types, expr_presence, datasets, contrast_types, default_display_cols, optional_display_cols, 
    default_stat_fields, optional_stat_fields) {
    
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
                    
                    # checkboxInput(ns("show_settings"), "Show settings", value=TRUE),
                    # conditionalPanel(
                        # sprintf("input['%s'] == 1", ns("show_settings")),

                    selectInput(ns("dataset"), "Dataset", choices=names(datasets), selected=names(datasets)[1]),
                    selectInput(ns("timepoint"), "Timepoint", choices=timepoints, selected="4d"),
                    selectInput(ns("contrast_type"), "Contrast type", choices=contrast_types, selected="Infection"),

                    checkboxInput(ns("do_fdr_filter"), "Do FDR and fold filtering", value=FALSE),
                    conditionalPanel(
                        sprintf("input['%s'] == 1", ns("do_fdr_filter")),
                        selectInput(ns("reg_type"), "Regulation", choices=c("all", "same", "contra"), selected="all"),
                        sliderInput(ns("fdr_cutoff_argamak"), "FDR cutoff cond. 1", value=0.1, step=0.01, min=0, max=1),
                        sliderInput(ns("fdr_cutoff_belinda"), "FDR cutoff cond. 2", value=0.1, step=0.01, min=0, max=1)
                    ),

                    checkboxInput(ns("trunc_long"), "Truncate long strings", value=TRUE),

                    checkboxInput(ns("do_annot_filter"), "Do annotation filtering", value=FALSE),
                    conditionalPanel(
                        sprintf("input['%s'] == 1", ns("do_annot_filter")),
                        selectInput(ns("annot_type"), "Annotation presence", choices=annot_types, selected=annot_types[1]),
                        selectInput(ns("arg_expr_pres"), "Cond. 1 expression presence", choices=expr_presence, selected="ALL"),
                        selectInput(ns("bel_expr_pres"), "Cond. 2 expression presence", choices=expr_presence, selected="ALL")
                    ),

                    selectInput(ns("table_add_shown_fields"), "Additional shown fields", choices=optional_display_cols, selected=default_display_cols, multiple=TRUE),
                    selectInput(ns("table_add_stat_fields"), "Additional stat fields", choices=optional_stat_fields, selected=default_stat_fields, multiple=TRUE),

                    actionButton(ns("button_show_align"), "Do Align"),
                    # textInput(ns("download_base_name"), "Download name", value="current_data"),
                    downloadButton(ns("download_current"), "Download selection"),
                    textOutput(ns("table_enrichment_status"))
                    # )
                ),
                fluidPage(
                    style = "flex-grow:1; resize:horizontal; overflow-x: scroll; overflow-y: hidden;",
                    fluidRow(
                        fluidPage(
                            h4("About"),
                            div("Filtering can be performed on FDR (infected-control), presence in transcriptome assemblies for peptides and protein expression"),
                            div("Do alignment by clicking row and press 'Do align'"),
                            div("Enrichment is performed for current filtering selection, with all (non-filtered) IDs as the universe")
                        ),
                        DT::dataTableOutput(ns("table"))
                    )
                )
            )
        )
    )
}

table_vis <- function(input, output, session, table_vars) {

    output$table <- DT::renderDataTable({
        show_table(
            table_vars$cached_filtered_table(),
            table_vars$stat_bases(),
            get_all_cols=FALSE,
            annot_cols=input$table_add_shown_fields,
            stat_cols=input$table_add_stat_fields
        )
    })
}

table_panel <- function(input, output, session, datasets, open_tab, sample_name="old_sample") {

    observeEvent(input$button_show_align, {
        open_tab("Alignment")
    })

    table_vars <- list()
    table_vars$stat_bases <- reactive({
        if (input$contrast_type == "Infection") 
            stat_base <- paste(c("Inf", "Ctl"), input$timepoint, sep="_")
        else if (input$contrast_type == "Variety") 
            stat_base <- paste(c("Arg", "Bel"), input$timepoint, sep="_")
        else 
            stop("Unknown contrast type: ", input$contrast_type)
        stat_base
    })
    
    filtered_table <- reactive({
        get_filter_table(
            datasets[[input$dataset]],
            table_vars$stat_bases(),
            fold_type=input$reg_type,
            fdr_cutoff_arg=input$fdr_cutoff_argamak,
            fdr_cutoff_bel=input$fdr_cutoff_belinda,
            annotation_presence=input$annot_type,
            argamak_expr=input$arg_expr_pres,
            belinda_expr=input$bel_expr_pres,
            do_fdr_filter=input$do_fdr_filter,
            do_string_truncate=input$trunc_long,
            include_sdf=TRUE,
            contrast_type=input$contrast_type
    )})
    
    get_settings <- function() {
        settings <- list()
        settings[["dataset"]] <- input$dataset
        settings[["do_fdr_filter"]] <- input$do_fdr_filter
        settings[["fdr_cutoff_argamak"]] <- input$fdr_cutoff_argamak
        settings[["fdr_cutoff_belinda"]] <- input$fdr_cutoff_argamak
        settings[["annot_type"]] <- input$annot_type
        settings[["stat_bases_1"]] <- table_vars$stat_bases()[1]
        settings[["stat_bases_2"]] <- table_vars$stat_bases()[2]
        settings_df <- do.call("rbind", settings) %>% data.frame()
        settings_df <- cbind(rownames(settings_df), settings_df)
        colnames(settings_df) <- c("Parameter", "Value")
        settings_df
    }
    
    current_settings <- reactive({
        get_settings()
    })

    output$download_current <- downloadHandler(
        filename = function() { 
            sprintf("%s.tsv", "current_data")
        },
        content = function(fname) {
            readr::write_tsv(filtered_table(), fname)
        }
    )
    
    table_vars$cached_full_table <- reactive({
        cbind(SummarizedExperiment::rowData(datasets[[input$dataset]]) %>% 
                  data.frame(), assay(datasets[[input$dataset]]) %>% data.frame())
    })
    
    table_vars$cached_filtered_table <- filtered_table
    table_vars$cached_sdf <- reactive({
        filtered_table() %>% dplyr::select(SummarizedExperiment::colData(datasets[[input$dataset]])[[sample_name]])
    })
    
    table_vars$dataset <- reactive({
        datasets[[input$dataset]]
    })
    
    table_vars$timepoint <- reactive({
        input$timepoint
    })

    table_vars$target_id <- reactive({
        target_row <- input$table_rows_selected[1]
        target_id <- filtered_table()[target_row, ]$ProteinID
        target_id
    })
    
    table_vars$contrast_type <- reactive({
        input$contrast_type
    })
    
    return(table_vars)
}

# Here - insert the three filtering aspects
get_filter_table <- function(dataset, stat_bases, fold_type="all", fdr_cutoff_arg=0.1, fdr_cutoff_bel=0.1,
                             annotation_presence="all", argamak_expr="ALL", belinda_expr="ALL", do_fdr_filter=TRUE,
                             do_string_truncate=FALSE, include_sdf=FALSE, contrast_type="Infection") {

    no_round_fields <- c(
        "EValue"
    )

    format_col <- function(col) {
        if (typeof(col) == "double" && !(col %in% no_round_fields)) {
            round(col, 5)
        }
        else if (typeof(col) == "character") {
            substr(col, 1, 20)
        }
        else {
            col
        }
    }

    fold_filter <- function(fold1, fold2, fold_type) {
        if (fold_type == "all") {
            TRUE
        }
        else if (fold_type == "same") {
            sign(fold1) == sign(fold2)
        }
        else if (fold_type == "contra") {
            sign(fold1) != sign(fold2)
        }
        else {
            stop("Unknown fold type: ", fold_type)
        }
    }

    filtered_df <- SummarizedExperiment::rowData(dataset) %>% data.frame()

    if (include_sdf) {
        filtered_df <- cbind(filtered_df, assay(dataset) %>% data.frame())
    }

    if (do_string_truncate) {
        filtered_df <- filtered_df %>%
            lapply(format_col) %>%
            data.frame()
    }

    if (do_fdr_filter) {
        filtered_df <- filtered_df %>%
            dplyr::filter(fold_filter(
                UQ(as.name(sprintf("%s.logFC", stat_bases[1]))),
                UQ(as.name(sprintf("%s.logFC", stat_bases[2]))), fold_type)) %>%
            dplyr::filter(
                UQ(as.name(sprintf("%s.adj.P.Val", stat_bases[1]))) < fdr_cutoff_arg &
                    UQ(as.name(sprintf("%s.adj.P.Val", stat_bases[2]))) < fdr_cutoff_bel)
    }

    if (annotation_presence != "all") {
        filtered_df <- filtered_df %>% dplyr::filter(.data$annot_type == annotation_presence)
    }

    if (argamak_expr != "ALL") {
        arg_filter_col <- sprintf("%s.presence", stat_bases[1])
        filtered_df <- filtered_df %>% dplyr::filter(UQ(as.name(arg_filter_col)) == argamak_expr)
    }

    if (belinda_expr != "ALL") {
        bel_filter_col <- sprintf("%s.presence", stat_bases[2])
        filtered_df <- filtered_df %>% dplyr::filter(UQ(as.name(bel_filter_col)) == belinda_expr)
    }
    
    filtered_df
}

show_table <- function(filtered_df, stat_bases, stat_cols=NULL, annot_cols=NULL, get_all_cols=FALSE, default_length=25) {

    if (!get_all_cols) {
        target_fields <- c(
            annot_cols,
            paste(stat_bases[1], stat_cols, sep="."),
            paste(stat_bases[2], stat_cols, sep=".")
        )
    }
    else {
        target_fields <- colnames(filtered_df)
    }
    
    select_filter_df <- filtered_df %>% dplyr::select(target_fields)
    select_filter_df %>%
        DT::datatable(selection='single', class="compact cell-border", options=list(pageLength=default_length)) %>%
        DT::formatStyle(
            c(sprintf("%s.logFC", stat_bases[1]), sprintf("%s.logFC", stat_bases[2])), 
            color = htmlwidgets::JS("value > 0 ? 'red': 'blue'")
        ) %>%
        DT::formatStyle(
            colnames(select_filter_df),
            fontSize = '80%'
        )
}

