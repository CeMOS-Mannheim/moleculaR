# Both tabs use separate session-owned Cardinal modules.
ui <- moleculaR:::.cardinalAppUi()
server <- function(input, output, session) {
    modules <- moleculaR:::.cardinalAppServer(input, output, session)
}
shiny::shinyApp(ui = ui, server = server)
