in_package <- F

if (!in_package) {
    source("module_enrichment.R")
    source("module_multivar_vis.R")
    source("module_about.R")
    source("module_table.R")
    source("module_single_feature.R")
    
    print(getwd())
    
    shiny::shinyOptions(
        query_proteins_arg_fp = "../../FusariumResponseInOatMethods_files/shiny_data/transcripts-shortid_fa_transdecoder-Arg-shortid_pep.fasta",
        query_proteins_bel_fp = "../../FusariumResponseInOatMethods_files/shiny_data/transcripts-shortid_fa_transdecoder-Bel-shortid_pep.fasta",
        search_fasta_fp = "../../FusariumResponseInOatMethods_files/shiny_data/search_protein.fasta",
        rds_obj_fp = "../../FusariumResponseInOatMethods_files/shiny_data/combined_flat_ses.rds"
    )
}

get_global <- function() {

    global <- list()
    
    # Libraries
    library(R6)
    library(shiny)
    library(tidyverse)
    library(SummarizedExperiment)
    library(DT)
    library(shinythemes)
    library(enrichplot)
    library(limma)
    library(clusterProfiler)
    library(enrichplot)
    library(msaR)
    library(org.At.tair.db)
    library(ggdendro)

    # Paths
    
    query_proteins_arg_fp <- shiny::getShinyOption("query_proteins_arg_fp")
    query_proteins_bel_fp <- shiny::getShinyOption("query_proteins_bel_fp")
    search_fasta_fp <- shiny::getShinyOption("search_fasta_fp")
    rds_obj_fp <- shiny::getShinyOption("rds_obj_fp")
    
    if (is.null(rds_obj_fp)) {
        message("No rds_obj_fp found, returning default settings")
        global$datasets <- "datasets"
        global$features <- "features"
        global$all_cols <- "all_cols"
        global$optional_display_cols <- "optional_display_cols"
        global$default_display_cols <- "default_display_cols"
        global$base_target <- "base_target"
        global$optional_stat_fields <- "optional_stat_fields"
        global$default_stat_fields <- "default_stat_fields"
        global$timepoints <- "timepoints"
        global$contrast_types <- "contrast_types"
        global$annot_types <- "annot_types"
        global$expr_presence <- "expr_presence"
        global$conditions <- "conditions"
        return(global)
    }
    
    datasets <- readRDS(rds_obj_fp)
    selected_dataset <- datasets[[1]]
    
    message("Loading query proteins")
    query_proteins_arg <- Biostrings::readAAStringSet(query_proteins_arg_fp)
    query_proteins_bel <- Biostrings::readAAStringSet(query_proteins_bel_fp)
    query_proteins <- c(query_proteins_arg, query_proteins_bel)
    
    load_string_set <- function(fp, name_delims=" ") {
        
        search_strings <- Biostrings::readAAStringSet(fp)
        search_names <- names(search_strings)
        for (name_delim in name_delims) {
            search_names <- limma::strsplit2(search_names, name_delim)[, 1]
        }
        
        names(search_strings) <- search_names
        search_strings
    }
    
    message("Load reference proteins, trimming out annotation part in both Arabidopsis and fungal")
    search_strings <- load_string_set(search_fasta_fp, name_delims=c(" ", "\\|"))
    message("Done loading!")
    
    features <- rowData(selected_dataset) %>% 
        data.frame() %>% 
        filter(!is.na(ProteinID)) %>% 
        dplyr::select(ProteinID) %>% 
        unlist() %>% 
        as.character() %>% 
        unname() %>% 
        sort() %>% 
        unique() %>% 
        strsplit(",") %>% 
        unlist() %>% 
        gsub("\\|.*", "", .) %>% 
        unique() %>% 
        sort()
    
    # Setup available filtering cols (preliminary solution)
    all_cols <- colnames(rowData(selected_dataset))
    optional_display_cols <- all_cols[!grepl("4h|1d|2d|4d", all_cols)]
    default_display_cols <- c("annot_type", "Protein", "ProteinID", "EValue", "Description")
    
    # Setup stat-fields (targeting the current timepoint)
    base_target <- "Arg_4d\\."
    optional_stat_fields <- gsub(base_target, "", all_cols[grepl(base_target, all_cols)])
    default_stat_fields <- c("AveExpr", "adj.P.Val", "logFC", "presence")
    
    # Additional menu parameter options
    timepoints <- c("4h", "1d", "2d", "4d")
    contrast_types <- c("Infection", "Variety")
    annot_types <- c("all", unique(rowData(datasets[[1]])$annot_type))
    expr_presence <- c("ALL", "BOTH", "NONE", "HIGHONLY", "LOWONLY")
    
    conditions <- colnames(colData(datasets[[1]]))
    
    global$datasets <- datasets
    global$features <- features
    global$all_cols <- all_cols
    global$optional_display_cols <- optional_display_cols
    global$default_display_cols <- default_display_cols
    global$base_target <- base_target
    global$optional_stat_fields <- optional_stat_fields
    global$default_stat_fields <- default_stat_fields
    global$timepoints <- timepoints
    global$contrast_types <- contrast_types
    global$annot_types <- annot_types
    global$expr_presence <- expr_presence
    global$conditions <- conditions
    
    global
}


