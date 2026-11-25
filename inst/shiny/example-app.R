ui <- moleculaR:::.cardinalAppUi(exampleOnly = TRUE)
server <- function(input, output, session) {
    modules <- moleculaR:::.cardinalAppServer(input, output, session, exampleOnly = TRUE)
}
shiny::shinyApp(ui = ui, server = server)
