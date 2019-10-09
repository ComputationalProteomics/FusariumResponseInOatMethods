about_panel_ui <- function(id) {
    ns <- NS(id)
    
    tabPanel(
        id,
        fluidPage(
            "Developed by Jakob Willforss at the Department of Immunotechnology, Lund University"
        )
    )
}