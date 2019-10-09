library(shiny)
library(tidyverse)

single_feature_panel_ui <- function(id, features, datasets_names, conditions) {
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
                    selectInput(ns("feature"), "Feature", choices=features, selected=features[1], width="300px"),
                    selectInput(ns("subid"), "SubID", choices=c(), selected=NULL, width="300px"),
                    
                    conditionalPanel(
                        sprintf("input['%s'] == 'Intensity'", ns("plot_tabs")),
                        selectInput(ns("datasetA"), "Dataset A", choices=datasets_names, selected=datasets_names[1]),
                        selectInput(ns("datasetB"), "Dataset B", choices=c("null", datasets_names), selected="null")
                    ),
                    
                    conditionalPanel(
                        sprintf("input['%s'] == 'Contrast'", ns("plot_tabs")),
                        selectInput(ns("contrast_dataset"), "Dataset", choices=datasets_names, selected=datasets_names[1]),
                        selectInput(ns("split_cond"), "Split cond.", choices=conditions, selected="comb_short", width="300px"),
                        checkboxInput(ns("only_4d"), "Only 4d", value=TRUE)
                    )
                ),
                fluidPage(
                    style = "flex-grow:1; resize:horizontal; overflow-x: hidden; overflow-y: hidden;",
                    tabsetPanel(
                        type="tabs",
                        id=ns("plot_tabs"),
                        tabPanel(
                            "Alignment",
                            textOutput(ns("mutation_status")),
                            tabPanel("MSAR", msaR::msaROutput(ns("alignment_msar")))
                        ),
                        tabPanel(
                            "Contrast",
                            plotOutput(ns("ContrastPlots"), height=800)
                        ),
                        tabPanel(
                            "Intensity",
                            plotOutput(ns("IntensityPlots"), height=800)
                        )
                    )
                )
            )
        )
    )
}

#' @importFrom SummarizedExperiment colData rowData
single_feature_panel <- function(input, output, session, table_vars, datasets, query_proteins, search_strings) {
    
    align <- reactiveVal(NULL)
    positions <- reactiveVal(NULL)

    output$alignment_msar <- msaR::renderMsaR({
        align(make_proteogenomic_alignment(
            table_vars$dataset(),
            c(query_proteins, search_strings),
            input$feature,
            input$dimension[1]
        ))
        positions(analyze_positions(align(), pattern1="^Ag", pattern2="^Be"))
        message("Identified positions: ", paste(positions(), collapse=", "))

        if (!is.null(align())) {
            msaR::msaR(align(), alignmentHeight=500, width=800, colorscheme="clustal2", overviewboxWidth="auto")
        }
    })

    output$mutation_status <- renderText({
        paste("Putative mutation sites:", paste(positions(), collapse=", "))
    })

    observeEvent(table_vars$target_id(), {
        if (length(table_vars$target_id()) != 0) {
            updateSelectInput(session, "feature", selected=table_vars$target_id())
        }
    })

    observeEvent(input$feature, {
        df <- rowData(table_vars$dataset()) %>% data.frame()
        entries <- df[grepl(input$feature, df$ProteinID), ]$External.IDs
        updateSelectInput(session, "subid", choices = entries)
    })

    output$IntensityPlots <- renderPlot({

        if (input$datasetA != input$datasetB && input$datasetB != "null") {

            do_intensity_plot(
                datasets[[input$datasetA]],
                feature=input$feature,
                datasetB=datasets[[input$datasetB]],
                assembly_id=input$subid
            )
        }
        else {
            do_intensity_plot(
                datasets[[input$datasetA]],
                input$feature,
                assembly_id=input$subid
            )
        }
    })

    output$ContrastPlots <- renderPlot({

        do_contrast_plot(
            datasets[[input$contrast_dataset]],
            feature=input$feature,
            assembly_id=input$subid,
            split_condition=input$split_cond,
            only_4d=input$only_4d
        )
    })
}

# Contrast panel

#' @import ggplot2
do_contrast_plot <- function(dataset, feature, assembly_id, split_condition, protein_id_name="ProteinID", sub_id_name="External.IDs", only_4d=FALSE, verbose=FALSE) {
    
    verbose <- TRUE
    if (verbose) {
        message("Visualizing feature: ", feature, " assembly_id: ", assembly_id, 
                " split_condition: ", split_condition, " protein_id_name: ", protein_id_name, 
                " sub_id_name: ", sub_id_name)
    }
    
    conds <- dataset %>% 
        colData() %>% 
        data.frame() %>% 
        dplyr::select(split_condition) %>% 
        unlist() %>% 
        unname()
    protein_col <- dataset %>% 
        rowData() %>% 
        data.frame() %>% 
        dplyr::filter(UQ(as.name(protein_id_name)) == feature)
    
    se_slice <- dataset[which(rowData(dataset)[[sub_id_name]] == assembly_id), ]
    
    if (only_4d) {
        se_slice <- se_slice[, SummarizedExperiment::colData(se_slice)$time == "4d"]
    }
    
    plot_df <- cbind(
        cond=SummarizedExperiment::colData(se_slice)[[split_condition]],
        assay(se_slice) %>% t() %>% data.frame()
    ) %>% 
        tidyr::gather("feature", "value", -.data$cond)
    
    ggplot(
        plot_df, 
        aes(x=.data$cond, y=.data$value, color=.data$feature)) + 
            geom_boxplot() + 
            geom_point(position=position_jitterdodge(jitter.width=0.05), size=3) + 
            ggtitle(assembly_id)
}

# Intensity panel

parse_dataset <- function(dataset, feature, feature_label, assembly_id=NULL, protein_id_name="ProteinID", sub_id_name="External.IDs") {
    
    parsed_ids <- limma::strsplit2(rowData(dataset)[[protein_id_name]], "\\|")[, 1]
    assembly_ids <- rowData(dataset)[[sub_id_name]]
    
    if (!is.null(assembly_id)) {
        matching_inds <- which(feature == parsed_ids & assembly_id == assembly_ids)
    }
    else {
        matching_inds <- which(feature == parsed_ids)
    }
    
    sub_assay <- assay(dataset)[matching_inds, , drop=FALSE] %>% data.frame()
    sub_assay$feature <- paste0(feature_label, seq_len(nrow(sub_assay)))
    sub_assay$group <- feature_label
    long_df <- sub_assay %>% 
        tidyr::gather("sample", "intensity", -"feature", -"group")
    long_df
}

parse_to_long_df <- function(feature, datasetA, datasetB, assembly_id=NULL) {
    
    longA <- parse_dataset(datasetA, feature, assembly_id=assembly_id, feature_label="A")
    if (!is.null(datasetB)) {
        longB <- parse_dataset(datasetB, feature, assembly_id=assembly_id, feature_label="B")
        long_combined <- rbind(longA, longB)
    }
    else {
        long_combined <- longA
    }
    
    long_combined    
}

#' @import ggplot2
do_intensity_plot <- function(datasetA, feature, assembly_id=NULL, datasetB=NULL) {
    
    long_combined <- parse_to_long_df(feature, datasetA, datasetB, assembly_id=assembly_id)
    ggplot(long_combined, aes(x=sample, y=.data$intensity, group=feature, color=feature)) + 
        ggtitle("Intensities") + 
        geom_point(na.rm=TRUE) + 
        geom_line(aes(linetype=.data$group))
}

# Alignment panel
#' @importFrom rlang .data
analyze_positions <- function(align, pattern1, pattern2, min_support=2) {

    get_site_string <- function(site_mat, pattern1, pattern2) {
        pat1_mat <- site_mat[grepl(pattern1, rownames(site_mat)), , drop=FALSE]
        pat2_mat <- site_mat[grepl(pattern2, rownames(site_mat)), , drop=FALSE]
        pat1_uniques <- apply(pat1_mat, 2, function(col) { unique(col[col != "-"]) })
        pat2_uniques <- apply(pat2_mat, 2, function(col) { unique(col[col != "-"]) })
        paste(pat1_uniques, pat2_uniques, sep="/")
    }
    
    is_pos_mutated <- function(test_col, pattern1, pattern2) {
        pat1_col <- test_col[grepl(pattern1, names(test_col))]
        pat2_col <- test_col[grepl(pattern2, names(test_col))]
        
        pat1_col_trimmed <- pat1_col[which(pat1_col != "-")]
        pat2_col_trimmed <- pat2_col[which(pat2_col != "-")]
        
        if (length(pat1_col_trimmed) > 0 && length(pat2_col_trimmed) > 0) {
            pat1_unique <- unique(pat1_col_trimmed)
            pat2_unique <- unique(pat2_col_trimmed)
            
            if (length(pat1_unique) == 1 && length(pat2_unique) == 1 && length(pat1_col_trimmed) > min_support && length(pat2_col_trimmed) > min_support) {
                pat1_unique != pat2_unique
            }
            else {
                FALSE
            }
        }
        else {
            FALSE
        }
    }
    
    align_mat <- as.matrix(align)
    mutation_contrast <- apply(align_mat, 2, is_pos_mutated, pattern1=pattern1, pattern2=pattern2)
    indices <- which(mutation_contrast)
    
    if (length(indices) > 0) {
        
        status_string <- get_site_string(align_mat[, indices, drop=FALSE], pattern1, pattern2)
        paste(indices, status_string, sep="-")
    }
    else {
        "NOTHING FOUND"
    }
}

make_proteogenomic_alignment <- function(dataset, search_sequences, target_feature, width, prompt_val=NULL, protein_id_col="ProteinID", pep_seq_col="Peptide.Sequence") {
    
    blast_hit_protein_ids <- limma::strsplit2(rowData(dataset)$ProteinID, "\\|")[, 1]
    matching_external_ids <- rowData(dataset)$External.IDs[blast_hit_protein_ids %in% target_feature]
    transcript_ids <- sort(unique(unlist(strsplit(matching_external_ids, ","))))

    if (pep_seq_col %in% colnames(rowData(dataset))) {
        rowData(dataset)$clean_peps <- get_clean_peptides(dataset, pep_seq_col)
        ms_peps <- rowData(dataset) %>%
            data.frame() %>%
            dplyr::filter(.data$ProteinID == target_feature) %>%
            dplyr::select(.data$clean_peps) %>%
            unlist() %>%
            Biostrings::AAStringSet()
    }
    else {
        ms_peps <- NULL
    }

    target_query_proteins <- search_sequences[transcript_ids]
    seqs <- Biostrings::AAStringSet(c(
        search_sequences[target_feature],
        target_query_proteins[sort(names(target_query_proteins))],
        ms_peps
    ))
    align <- DECIPHER::AlignSeqs(seqs)
    align
}

get_clean_peptides <- function(dataset, peptide_col) {
    if (!(peptide_col %in% colnames(rowData(dataset)))) {
        NULL
    }
    else {
        pep_seqs <- rowData(dataset)[[peptide_col]]
        pep_seqs_splits <- limma::strsplit2(pep_seqs, ",| ")
        clean_pep_seqs <- pep_seqs_splits[, 1]
        clean_pep_seqs        
    }
}



