dirpath <- file.path("~", "founder_diet_study")
dirpath <- file.path(dirpath, "HarmonizedData")
traitSignal <- readRDS(file.path(dirpath, "liverSignal.rds"))
traitStats <- readRDS(file.path(dirpath, "liverStats.rds"))

################################################################

title <- "Test Shiny Module"

ui <- function() {
  
  shiny::fluidPage(
    shiny::titlePanel(title),
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        # Datasets and Traits.
        shiny::fluidRow(
          shiny::column(4, shiny::uiOutput("dataset")),
          shiny::column(4, shiny::selectInput("sex", "Sex:",
            c("Both Sexes", "Female", "Male", "Sex Contrast"))),
          shiny::column(4, shiny::uiOutput("term"))),
        
        shiny::fluidRow(
          shiny::column(4, shiny::sliderInput("p.value", "p.value:",
                                              0, 1, 0.05, 0.05)),
          shiny::column(4, shiny::sliderInput("power", "power:",
                                              6, 30, 6, 6)),
          shiny::column(4, shiny::sliderInput("size", "size:",
                                              2, 14, 4, 2))),
        
        shiny::fluidRow(
          shiny::column(6, shiny::selectInput("varyParam", "Vary:",
                                              c("size","power"))),
          shiny::column(6, shiny::checkboxInput("isDendro", "Dendrogram?")))
        ),
      
      shiny::mainPanel(
        shiny::uiOutput("intro"),
        shiny::uiOutput("dendro"),
        shiny::uiOutput("topPower"),
        )
    ))
}

server <- function(input, output, session) {
  # Related Datasets.
  
  output$dataset <- renderUI({
    datasets <- unique(traitStats$dataset)
    shiny::selectInput("dataset", "Datasets:", datasets, "LivMet",
                       multiple = TRUE)
  })
  output$term <- renderUI({
    terms <- unique(traitStats$term)
    shiny::selectInput("term", "Term:", terms, terms[1])
  })
  
  traitSignalInput <- reactive({
    shiny::req(input$dataset)
    
    dplyr::filter(traitSignal, dataset %in% input$dataset)
  })
  traitStatsInput <- reactive({
    shiny::req(input$dataset)
    
    dplyr::filter(traitStats, dataset %in% input$dataset)
  })
  traitContr <- reactive({
    shiny::req(traitSignalInput(), traitStatsInput(), input$term)
    
    foundr::conditionContrasts(traitSignalInput(), traitStatsInput(),
                                 termname = input$term)
  })
  traitContrPval <- reactive({
    shiny::req(traitContr(), input$p.value)
    
    dplyr::filter(shiny::req(traitContr()), .data$p.value <= input$p.value)
  })
  
  output$intro <- renderUI({
    shiny::renderText("intro", {
      paste("Guideline is to have power of 6 and size of 4 for unsigned modules.",
            "Power larger than 30 is too high.",
            "Topology curves are not necessarily helpful.",
            "This tool lets one explore the topology a bit.")
    })
  })

  powers <- seq(from = 6, to = 30, by = 6)
#  powers <- seq(from = 10, to = 200, by = 10)
  sizes <- seq(from = 2, to = 14, by = 2)
  output$topPower <- renderUI({
    shiny::req(traitContrPval(), input$sex)
    
    top <- modulr::wgcna_topology(
      dplyr::filter(
        dplyr::select(traitContrPval(), -p.value),
        .data$sex == input$sex),
      power = powers)
    
    shiny::tagList(
      shiny::h3("WGCNA Topology"),
      shiny::renderPlot(print(modulr::ggplot_wgcna_topology(top)))
    )
  })
  
  dendro <- shiny::reactive({
    shiny::req(traitContrPval(), input$varyParam, input$isDendro)
    
    mod <- list()
    # More than 13 causes crash.
    shiny::withProgress(
      message = paste("Dendrogram", 'calculations in progress'),
      detail = 'This may take a while...',
      value = 0,
      {
        switch(
          input$varyParam,
          size = {
            shiny::req(input$power)
            
            for(size in sizes) {
              mod[[paste0("size=", size)]] <- tryCatch(
                modulr::listof_wgcnaModules(
                  dplyr::select(traitContrPval(), -p.value),
                  params = list(power = input$power, minSize = size)),
                error = function(e) NULL)
              shiny::incProgress(1/length(sizes))
            }},
          power = {
            shiny::req(input$size)
            
            for(beta in powers) {
              mod[[paste0("beta=", beta)]] <- tryCatch(
                modulr::listof_wgcnaModules(
                  dplyr::select(traitContrPval(), -p.value),
                  params = list(power = beta, minSize = input$size)),
                error = function(e) NULL)
              shiny::incProgress(1/length(powers))
            }})
      })
    
    if(!length(mod))
      return(NULL)

    purrr::transpose(mod)
  })
    
  output$dendro <- renderUI({
    shiny::req(dendro())
    
    response <- switch(
      shiny::req(input$varyParam),
      size = paste0("size=", shiny::req(input$size)),
      power = paste0("beta=", shiny::req(input$power)))
    
    p <- foundr::ggplot_listof_wgcnaModules(dendro()[[input$sex]], response)

    shiny::tagList(  
      shiny::h3("WGCNA Dendrograms"),
      shiny::renderPlot(print(p)))
  })
}

shiny::shinyApp(ui = ui, server = server)
