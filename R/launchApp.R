#' launches the shinyAppDemo app
#'
#' @export launchApp
#' @return shiny application object
#' @param query_proteins_arg_fp TODO
#' @param query_proteins_bel_fp TODO
#' @param search_fasta_fp TODO
#' @param rds_obj_fp TODO
#' @import shiny SummarizedExperiment dplyr rlang
launchApp <- function(
    query_proteins_arg_fp = "../../FusariumResponseInOatMethods_files/shiny_data/transcripts-shortid_fa_transdecoder-Arg-shortid_pep.fasta",
    query_proteins_bel_fp = "../../FusariumResponseInOatMethods_files/shiny_data/transcripts-shortid_fa_transdecoder-Bel-shortid_pep.fasta",
    search_fasta_fp = "../../FusariumResponseInOatMethods_files/shiny_data/search_protein.fasta",
    rds_obj_fp = "../../FusariumResponseInOatMethods_files/shiny_data/combined_flat_ses.rds") {
    
    shiny::shinyOptions(
        query_proteins_arg_fp = query_proteins_arg_fp,
        query_proteins_bel_fp = query_proteins_bel_fp,
        search_fasta_fp = search_fasta_fp,
        rds_obj_fp = rds_obj_fp
    )
    
    if (is.null(any(c(query_proteins_arg_fp, query_proteins_bel_fp, search_fasta_fp, rds_obj_fp)))) {
        stop(sprintf(
            "All four input files are required, obtained: \nquery_proteins_arg_fp: %s\nquery_proteins_bel_fp: %s\nsearch_fasta_fp: %s\nrds_obj_fp: %s",
            query_proteins_arg_fp, query_proteins_bel_fp, search_fasta_fp, rds_obj_fp
        ))
    }
    shinyApp(ui = ui, server = server)
}

