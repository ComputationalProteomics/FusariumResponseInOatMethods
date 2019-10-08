ui <- shinyUI(navbarPage(
    theme = shinytheme("flatly"),
    "OmicsLoupe",
    id="navbar",
    table_panel_ui("Table", timepoints, annot_types, expr_presence),
    pca_panel_ui("MultivarVis"),
    single_feature_panel_ui("SingleFeature", features, names(datasets)),
    enrichment_panel_ui("Enrichment"),
    about_panel_ui("About")
))
