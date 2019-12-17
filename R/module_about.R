about_panel_ui <- function(id) {
    ns <- NS(id)
    
    tabPanel(
        id,
        fluidPage(
            h3("About"),
            "Developed by Jakob Willforss at the Department of Immunotechnology, Lund University",
            h3("Help"),
            tabsetPanel(
                type="tabs",
                tabPanel(
                    "Table",
                    p("Displays a table illustrating the currently targeted set of data."),
                    p("Options: "),
                    uiOutput(ns("about_table"))
                ),
                tabPanel(
                    "MultivarVis",
                    p("Multivariate explorations of dataset using an interactive principal component analysis, 
                      hierarchical clustering or histogram of characteristics of interest. For the former two, subsets of the data can be 
                      studied, for example to investigate the dataset without outliers or to inspect certain subsets of the data"),
                    p("Options: "),
                    uiOutput(ns("about_multivar"))
                ),
                tabPanel(
                    "SingleFeature",
                    p("Univariate explorations of dataset, either on sequence level or expression level. The sequence alignment aligns
                      all assembly sequences homology matching to the target TAIR ID. The expression level allows exploration of intensity
                      levels across different sample conditions and inspects peptides matching to unique combinations of oat transcript IDs."),
                    p("Options: "),
                    uiOutput(ns("about_single"))
                ),
                tabPanel(
                    "Enrichment",
                    p("Perform enrichment of filtered subsets of the data. Either Gene Ontology Enrichment where the current subset
                      is compared to the overall dataset to see if certain GO-terms are overrepresented, or Gene Set Enrichment where ranked
                      fold changes are used for enrichment."),
                    p("Options: "),
                    uiOutput(ns("about_enrich"))
                )
            )
        )
    )
}

make_bullets <- function(text_vect) {
    sprintf("<ul><li>%s</li></ul>", paste(text_vect, collapse="</li><li>"))
}

about_panel <- function(input, output, session) { 
    
    ns <- session$ns
    
    output$about_table <- renderUI(
        HTML(
            make_bullets(c(
                "<b>Dataset</b> What dataset to inspect. The data is analyzed on peptide and protein level. Also, two batches in the dataset are analyzed both jointly and separately.", 
                "<b>Timepoint</b> Illustrate statistical comparisons for samples taken 4 hours, 1 day, 2 days and 4 days after infection",
                "<b>Contrast type</b> Inspect either contrasts between infected and control samples within varieties, or cross-variety comparisons",
                "<b>Do FDR and fold filtering</b> Options for FDR filtering the data<ul>",
                "<b>Regulation</b> Allows filtering out only features with same regulation or contra regulation",
                "<b>FDR cutoff cond. 1</b> FDR cutoff for the first selected statistical contrast",
                "<b>FDR cutoff cond. 2</b> FDR cutoff for the second selected statistical contrast</ul>",
                "<b>Truncate long strings</b> Truncate long strings in the table to prevent them overflowing the cells",
                "<b>Do annotation filtering</b> Filtering based on Argamak / Belinda presence in annotation and expression<ul>",
                "<b>Annotation presence</b> Filter for peptides matching to transcripts only from certain assemblies",
                "<b>Cond. 1 expression presence</b> Filter for expression presence in the first condition",
                "<b>Cond. 2 expression presence</b> Filter for expression presence in the second condition</ul>",
                "<b>Additional shown fields</b> Select columns of data to display (non-statistics)",
                "<b>Additional stat fields</b> Select columns of data to display (statistics, will display for both conditions)",
                "<b>[Button] Do Align</b> After marking a feature in the table, directly perform alignment for its reference transcripts",
                "<b>Download name</b> Base name when downloading the data",
                "<b>[Button] Download selection</b> Perform download of currently selected table data"
            ))
        )
    )
    
    output$about_multivar <- renderUI(
        HTML(
            make_bullets(c(
                "<b>Visualization</b> Type of illustration (PCA, hierarchical clustering or histogram)", 
                "<b>Condition</b> What sample condition to illustrate as main condition (used for coloring)", 
                "<b>Custom names</b> Assign custom labels to samples",
                "<b>Omit samples</b> Specify samples to omit (before recalculating overview plot)",
                "<b>Filter on type</b> Select sample condition to subset samples on",
                "<b>Levels to inspect</b> Related to previous condition - mark what levels to display",
                "<b>Show PCA settings</b> Toggle PCA settings <ul>",
                    "<b>Shape cond</b> Condition for shape of dots",
                    "<b>Variance filter</b> Remove low-variance features in the calculations",
                    "<b>Component 1</b> What component is PC1",
                    "<b>Component 2</b> What component is PC2",
                    "<b>Pairplot display</b> Display a grid of tiny PCA plots to screen for high-level patterns",
                    "<b>Pairplot count</b> Grid size if doing pairplot display</ul>",
                "<b>Show Scree</b> Show loadings of different principal components",
                "<b>Label size</b> Size of labels",
                "<b>Point size</b> Size of points"
            ))
        )
    )
    
    output$about_single <- renderUI(
        HTML(
            make_bullets(c(
                "<b>Feature</b> What TAIR homology ID to show", 
                "<b>Align only SubID IDs</b> Limits the aligned assembly sequences to those specifically matched by the peptide",
                "<b>SubID</b> For alignment: The the combination of assembly IDs matched by one peptide for which to perform sub alignment (if previous option is specified). For contrasts, what target transcript to inspect",
                "<b>Dataset</b> The target dataset to inspect (same as in Table tab)",
                "<b>Split cond.</b> Across over what sample condition to split the dataset (one condition level per x position)",
                "<b>Only 4d</b>Inspect only samples 4 days after infection (the main comparison inspected in the authors study)"
            ))
        )
    )
    
    output$about_enrich <- renderUI(
        HTML(
            make_bullets(c(
                "<b>[Button] Perform enrichment</b> Calculate enrichment for currently filtered part of dataset (in Table tab)", 
                "<b>Enrichment type</b> Either Gene Ontology Enrichment (GOE) or Gene Set Enrichment (GSEA)", 
                "<b>Ontology</b> The target ontology for which to calculate enrichment",
                "<b>Enrichment plot</b> The target plot type to illustrate the enrichment",
                "<b>Use BH FDR</b> For the enrichment, use FDR adjusted cutoff",
                "<b>FDR Cutoff</b> Cutoff for enrichment (FDR or p)",
                "<b>Displayed plot</b> Show plots for both comparisons, or only either",
                "<b>Plot cols</b> How many columns to show plots in",
                "<b>Plot height</b> Height of plots",
                "<b>Display terms</b> Number of enriched terms to illustrate"
            ))
        )
    )
}