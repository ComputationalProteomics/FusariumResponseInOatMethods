ui <- function() {
    shinyUI({
        message("Calling UI")
        global <- get_global()
        
        navbarPage(
            theme = shinythemes::shinytheme("flatly"),
            "OatOmics",
            id="navbar",
            table_panel_ui(
                "Table", 
                timepoints=global$timepoints, 
                annot_types=global$annot_types, 
                expr_presence=global$expr_presence, 
                datasets=global$datasets, 
                contrast_types=global$contrast_types,
                default_display_cols=global$default_display_cols, 
                optional_display_cols=global$optional_display_cols,
                default_stat_fields=global$default_stat_fields,
                optional_stat_fields=global$optional_stat_fields
            ),
            multivarvis_panel_ui(
                "MultivarVis",
                datasets=global$datasets,
                conditions=global$conditions
            ),
            single_feature_panel_ui(
                "SingleFeature", 
                features=global$features, 
                datasets_names=names(global$datasets),
                conditions=global$conditions
            ),
            enrichment_panel_ui("Enrichment"),
            about_panel_ui("About/Help")
        )
    })
}
