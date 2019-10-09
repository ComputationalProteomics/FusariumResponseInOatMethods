#' launches the shinyAppDemo app
#'
#' @export launchApp
#' @return shiny application object
#' @param query_proteins_arg_fp TODO
#' @param query_proteins_bel_fp TODO
#' @param search_fasta_fp TODO
#' @param rds_obj_fp TODO
#' @import shiny SummarizedExperiment dplyr rlang
launchApp <- function(query_proteins_arg_fp, query_proteins_bel_fp, search_fasta_fp, rds_obj_fp) {
    
    if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))
    
    # if (is.null(any(c(query_proteins_arg_fp, query_proteins_bel_fp, search_fasta_fp, rds_obj_fp)))) {
    #     stop(sprintf(
    #         "All four input files are required, obtained: \nquery_proteins_arg_fp: %s\nquery_proteins_bel_fp: %s\nsearch_fasta_fp: %s\nrds_obj_fp: %s",
    #         query_proteins_arg_fp, query_proteins_bel_fp, search_fasta_fp, rds_obj_fp
    #     ))
    # }
    # 
    # shiny::shinyOptions(in_path=in_path)
    shinyApp(ui = ui, server = server)
}

