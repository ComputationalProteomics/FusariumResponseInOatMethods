get_global <- function() {

    global <- list()
    
    query_proteins_arg_fp <- shiny::getShinyOption("query_proteins_arg_fp", NULL)
    query_proteins_bel_fp <- shiny::getShinyOption("query_proteins_bel_fp", NULL)
    search_fasta_fp <- shiny::getShinyOption("search_fasta_fp", NULL)
    rds_obj_fp <- shiny::getShinyOption("rds_obj_fp", NULL)
    
    message("Loaded: ", rds_obj_fp)
    
    if (is.null(rds_obj_fp)) {
        
        se <- SummarizedExperiment::SummarizedExperiment(
            assays=list(a=matrix(c(1,2,3,4), ncol=2)),
            colData=data.frame(a=c(1,2), b=c(2,3))
        )
        
        message("No rds_obj_fp found, returning default settings x")
        global$all_cols <- "all_cols"
        global$annot_types <- "annot_types"
        global$base_target <- "base_target"
        global$conditions <- "conditions"
        global$contrast_types <- "contrast_types"
        global$datasets <- list(a=se, b=se)
        global$default_display_cols <- "default_display_cols"
        global$default_stat_fields <- "default_stat_fields"
        global$expr_presence <- "expr_presence"
        global$features <- "features"
        global$optional_display_cols <- "optional_display_cols"
        global$optional_stat_fields <- "optional_stat_fields"
        global$query_proteins <- "query_proteins"
        global$search_strings <- "search_strings"
        global$timepoints <- "timepoints"
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
        dplyr::filter(!is.na(.data$ProteinID)) %>% 
        dplyr::select(.data$ProteinID) %>% 
        unlist() %>% 
        as.character() %>% 
        unname() %>% 
        sort() %>% 
        unique() %>% 
        strsplit(",") %>% 
        unlist() %>% 
        gsub("\\|.*", "", .data$.) %>% 
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
    
    conditions <- colnames(SummarizedExperiment::colData(datasets[[1]]))
    
    global$all_cols <- all_cols
    global$annot_types <- annot_types
    global$base_target <- base_target
    global$conditions <- conditions
    global$contrast_types <- contrast_types
    global$datasets <- datasets
    global$default_display_cols <- default_display_cols
    global$default_stat_fields <- default_stat_fields
    global$expr_presence <- expr_presence
    global$features <- features
    global$optional_display_cols <- optional_display_cols
    global$optional_stat_fields <- optional_stat_fields
    global$query_proteins <- query_proteins
    global$search_strings <- search_strings
    global$timepoints <- timepoints
    
    global
}


