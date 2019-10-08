#' launches the shinyAppDemo app
#'
#' @export launchApp
#' @return shiny application object
#' @example \dontrun {launchApp()}
#' @import shiny

library(shiny)

# wrapper for shiny::shinyApp()
launchApp <- function(in_path="dummy_path") {
    
    shiny::shinyOptions(in_path=in_path)
    
    shinyApp(ui = ui, server = server)
}
