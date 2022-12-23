library(shiny)
library(shinyWidgets)
# library(shinythemes)
# library(shinyjs)

# define variables in global environment
helpTxt1 <- paste0("The moleculaR R package provides a computational framework that introduces probabilistic 
                  mapping and point-for-point statistical testing of metabolites in tissue via Mass spectrometry 
                  imaging. It enables collective projections of metabolites and consequently spatially-resolved 
                  investigation of ion milieus, lipid pathways or user-defined biomolecular ensembles within the same image.")

helpTxt2 <- paste0("This companion shiny-app lets the user upload her own centroided imzML data and apply spatial probabilistic 
                  mapping through Molecular probabilistic Maps (MPMs) and Collective Projection probabilistic Maps (CPPMs). 
                  The first step would be to upload you own data. To do this click on `imzML & ibd` File and navigate to the 
                  directory containing the files you want to upload. Then select both imzML- and ibd-Files and upload them. 
                  Via `Spectrum .tsv File` uPload a continuous spectrum in `tsv` format which would represent either a single 
                  random pixel or a mean spectrum of your imaging dataset. Once this is done, please start the FWHM estimation
                  by presseing on 'Initialize'. For a given m/z value of interest, you can generate the corresponding MPM by 
                  providing that m/z value in the corresponding text box. To generate a collective projection probabilstic map (CPPM)
                  of a custom list of m/z values, please paste these values into the corresponding text box, comma separated.")                 
            



####  Frontend ####
ui <- navbarPage(p("moleculaR: Spatial Probabilistic Mapping of Metabolites in Mass Spectrometry Imaging", 
                                 HTML('&emsp;'), 
                                 HTML('&emsp;')),  
                        selected="Main",
                        
                        #theme = shinythemes::shinytheme("flatly"),
                        
                        tabPanel("Help",
                                    
                                 # tags$figure(
                                 #   align = "center",
                                 #   tags$img(src = "moleculaR-icon1.svg", 
                                 #            width = 500)
                                 # )
                                 includeMarkdown("./www/README.md")
                                #,
                                 #htmlOutput("helpTxt")
                        ),
                        
                        tabPanel("Main",
                                        tags$head(
                                          tags$style(HTML("hr {border-top: 1px solid #000000;}"))
                                          ),
                                        pageWithSidebar(
                                          headerPanel(""),
                                          
                                          sidebarPanel(
                                            wellPanel(
                                              fluidRow(
                                                column(12,
                                                              HTML(paste("<b>","Upload Files","</b>")),
                                                              fileInput("imzmlFile", "", 
                                                                               multiple = TRUE, 
                                                                               buttonLabel = "imzML & ibd Files",
                                                                               accept = c(".imzML",".ibd")),
                                                              fileInput("spectrFile", 
                                                                               NULL, 
                                                                               multiple = FALSE, 
                                                                               buttonLabel = "Spectrum .tsv File",
                                                                               accept = c(".tsv")),
                                                              actionButton(inputId = "loadData", 
                                                                                  label = "Load and Initialize",
                                                                                  style='padding:6px; font-size:80%')
                                                              
                                                ),
                                              )
                                            ),
                                            #hr(),
                                            
                                            # wellPanel(
                                            #   fluidRow(
                                            #   column(12,
                                            #                 HTML(paste0("<b>","Peaks FWHM Estimation","</b>")),
                                            #                 HTML(paste0("<b>"))
                                            #                 )),
                                            #   fluidRow(
                                            #     column(12, 
                                            #                   actionButton(inputId = "go_load", 
                                            #                                       label = "Initialize",
                                            #                                       style='padding:6px; font-size:80%')
                                            #                   )
                                            #   )),
                                         
                                        
                                            #hr(),
                                            
                                            #shinyjs::useShinyjs(),
                                            wellPanel(
                                              fluidRow(
                                                column(12,
                                                       HTML(paste0("<b>","Molecular Probabilitic Mapping","</b>"))),
                                                column(12,
                                                       numericInput("mz", label = "", value = 788.5447),
                                                       shinyWidgets::prettyCheckbox(inputId = "adjustmz",
                                                                                    label = "Find Closest Detectable m/z",
                                                                                    value = FALSE, 
                                                                                    icon = icon("check"), 
                                                                                    animation = "pulse"),
                                                       shinyWidgets::prettyCheckbox(inputId = "detailedPlotMPM",
                                                                                    label = "Show Detailed Plot", 
                                                                                    value = FALSE, 
                                                                                    icon = icon("check"), 
                                                                                    animation = "pulse"),
                                                              actionButton(inputId = "renderMPM", 
                                                                                  label = "Generate Plot",
                                                                                  style='padding:6px; font-size:80%'),
                                                ))
                                            ) ,
                                            
                                            
                                            wellPanel(
                                              fluidRow(
                                                column(12,
                                                              HTML(paste0("<b>","Collective Projection Probabilitic Mapping","</b>"))),
                                                column(12,
                                                              textInput("mzStr", 
                                                                        label = "Input comma separated m/z values", 
                                                                        placeholder  = "321.4321,123.1234,..."),
                                                       shinyWidgets::prettyCheckbox(inputId = "detailedPlotCPPM",
                                                                                    label = "Show Detailed Plot", 
                                                                                    value = FALSE, 
                                                                                    icon = icon("check"), 
                                                                                    animation = "pulse"),
                                                              
                                                              actionButton(inputId = "renderCPPM", 
                                                                                  label = "Generate Plot",
                                                                                  style='padding:6px; font-size:80%'),
                                                ))
                                            ),

                                            hr(),
                                            
                                            fluidRow(
                                              column(12,
                                                            actionButton(inputId = "exit", 
                                                                         label = "Exit App",
                                                                         style='padding:6px; font-size:80%', 
                                                                         width = "50%")
                                              )),
                                            width = 3),
                                          
                                          mainPanel(
                                            verbatimTextOutput("message"),
                                            plotOutput("imgs", width = 900, height=900) 
                                            
                                            )
                                          )),
                 tabPanel("Example",
                          tags$head(
                            tags$style(HTML("hr {border-top: 1px solid #000000;}"))
                          ),
                          pageWithSidebar(
                            headerPanel(""),
                            
                            sidebarPanel(
                              wellPanel(
                                fluidRow(
                                  column(12,
                                         HTML(paste("<b>","Load Example","</b>")),
                                         
                                         actionButton(inputId = "loadExample", 
                                                      label = "Load and Initialize",
                                                      style='padding:6px; font-size:80%')
                                         
                                  ),
                                )
                              ),
                              
          
                              wellPanel(
                                fluidRow(
                                  column(12,
                                         HTML(paste0("<b>","Molecular Probabilitic Mapping","</b>"))),
                                  column(12,
                                         numericInput("mzExample", label = "", value = 788.5447),
                                         shinyWidgets::prettyCheckbox(inputId = "adjustmzExample",
                                                                      label = "Find Closest Detectable m/z",
                                                                      value = TRUE, 
                                                                      icon = icon("check"), 
                                                                      animation = "pulse"),
                                         shinyWidgets::prettyCheckbox(inputId = "detailedPlotMPMExample",
                                                                      label = "Show Detailed Plot", 
                                                                      value = FALSE, 
                                                                      icon = icon("check"), 
                                                                      animation = "pulse"),
                                         actionButton(inputId = "renderMPMExample", 
                                                      label = "Generate Plot",
                                                      style='padding:6px; font-size:80%'),
                                  ))
                              ) ,
                              
                              
                              wellPanel(
                                fluidRow(
                                  column(12,
                                         HTML(paste0("<b>","Collective Projection Probabilitic Mapping","</b>"))),
                                  column(12,
                                         textInput("mzStrExample", 
                                                   label = "Input comma separated m/z values", 
                                                   placeholder  = "321.4321,123.1234,..."),
                                         shinyWidgets::prettyCheckbox(inputId = "adjustmzExample",
                                                                      label = "Find Closest Detectable m/z",
                                                                      value = TRUE, 
                                                                      icon = icon("check"), 
                                                                      animation = "pulse"),
                                         shinyWidgets::prettyCheckbox(inputId = "detailedPlotCPPMExample",
                                                                      label = "Show Detailed Plot", 
                                                                      value = FALSE, 
                                                                      icon = icon("check"), 
                                                                      animation = "pulse"),
                                         
                                         actionButton(inputId = "renderCPPMExample", 
                                                      label = "Generate Plot",
                                                      style='padding:6px; font-size:80%'),
                                  ))
                              ),
                              
                              hr(),
                              
                              fluidRow(
                                column(12,
                                       actionButton(inputId = "exitExample", 
                                                    label = "Exit App",
                                                    style='padding:6px; font-size:80%', 
                                                    width = "50%")
                                )),
                              width = 3),
                            
                            mainPanel(
                              verbatimTextOutput("messageExample"),
                              plotOutput("imgsExample", width = 900, height=900) 
                              
                            )
                          ))
                        
                        
                        )



####Backend ####

server <- function(input, output, session) {
  
  # increase maximum file uploade size to 48 GB
  options(shiny.maxRequestSize=48000*1024^2)
  
  
  # observeEvent for loading data
  observeEvent(input$loadData,{


    # checks
    output$message <- renderText({
      validate(
        need(length(input$spectrFile$datapath) != 0, "A spectrum has not been uploaded.")
      )
    })
    
    validate(
      need(length(input$spectrFile$datapath) != 0, "A spectrum has not been uploaded.")
    )
    
    output$message <- renderText({
      validate(
        need(grepl(".tsv", input$spectrFile$datapath), "The .tsv file (spectrum) has not been uploaded.")
      )
    })
    validate(
      need(grepl(".tsv", input$spectrFile$datapath), "The .tsv file (spectrum) has not been uploaded.")
    )
    

    output$message <- renderText({
      validate(
        need(any(grepl(".imzML", input$imzmlFile$datapath)), "The .imzMl file has not been uploaded.")
      )
    })
    validate(
      need(any(grepl(".imzML", input$imzmlFile$datapath)), "The .imzMl file has not been uploaded.")
    )

    output$message <- renderText({
      validate(
        need(any(grepl(".ibd", input$imzmlFile$datapath)), "The .ibd file has not been uploaded.")
      )
    })
    validate(
      need(any(grepl(".ibd", input$imzmlFile$datapath)), "The .ibd file has not been uploaded.")
    )

    # unify ibd and imzml names
    idx <- grep(".imzML", input$imzmlFile$datapath)
    oldName <- input$imzmlFile$datapath[idx]
    newName <- file.path(dirname(input$imzmlFile$datapath[idx]),"msiData.imzML")
    file.rename(from = oldName, to = newName)
    #input$imzmlFile$datapath[idx] <- newName
    
    idx <- grep(".ibd", input$imzmlFile$datapath)
    oldName <- input$imzmlFile$datapath[idx]
    newName <- file.path(dirname(input$imzmlFile$datapath[idx]),"msiData.ibd")
    file.rename(from = oldName, to = newName)
    #input$imzmlFile$datapath[idx] <- newName
  
    
    # start loading
    withProgress(
                message="please wait",
                detail="Loading and Processing Data...",
                value=0.2,{ 
                  
                  n <- 6
                  
                  # read the single spectrum
                  incProgress(1/n, detail = paste("Reading single spectrum .. "))
                  s          <<- moleculaR::readSingleSpect(input$spectrFile$datapath) 
                  
                  
                  # creat fwhm object
                  incProgress(1/n, detail = paste("Estimating FWHM .. "))
                  fwhmObj           <<- moleculaR::estimateFwhm(s = s) 
                  
                  # read the imzml data
                  incProgress(1/n, detail = paste("Reading imzML data .. "))
                  #imzmlPath <- input$imzmlFile$datapath[grep(".imzML", input$imzmlFile$datapath)]
                  
                  f <- list.files(path = dirname(input$imzmlFile$datapath), full.names = TRUE)
                  imzmlPath <- grep(pattern = ".imzML", x = f, value = TRUE)[1]
                  msData            <<- moleculaR::readCentrData(path = imzmlPath)
                  

                  # bin peaks
                  incProgress(1/n, detail = paste("Binning peaks .. "))
                  msData            <<-  MALDIquant::binPeaks(msData, 
                                                             tolerance = moleculaR::getFwhm(fwhmObj, 400)/400, #focusing on lipids 
                                                             method = "relaxed")
                  
                  # create sparse matrix representation
                  incProgress(1/n, detail = paste("Creating sparse matrix representation .. "))
                  spData <<- moleculaR::createSparseMat(x = msData)
                  
                  # create a spatial window
                  incProgress(1/n, detail = paste("Creating a spatial window .. "))
                  
                  output$imgs <- renderPlot({
                    par(mfrow = c(2, 1))
                    spwin     <<- moleculaR::createSpatialWindow(pixelCoords = MALDIquant::coordinates(msData), 
                                                                 clean = TRUE,
                                                                 plot = TRUE)
                    
                    plot(fwhmObj)
                  })
                  
                  
                  })
                
    
  })
  
  # observeEvent for MPM
  observeEvent(input$renderMPM,{
    
    
    # checks
    if(!exists("fwhmObj")){
      
      output$message <- renderText({
        "No FWHM data is available. Please upload and load data first."
      })
      
      output$imgs <- renderPlot({
        plot.new()
      })
      
      validate(
        need(exists("fwhmObj"), "No FWHM data is available. Please upload and load data first.")
      )
      
    }
    
    if(!exists("msData")){
      
      output$message <- renderText({
        "MSI data is not loaded. Please upload and load data first."
      })
      
      output$imgs <- renderPlot({
        plot.new()
      })
      
      validate(
        need(exists("msData"), "MSI data is not loaded. Please upload and load data first.")
      )
      
    }
    
    
    
    if(!exists("spData")){
      output$message <- renderText({
        "sparse data represention is not available. Please upload and load data first."
      })
      output$imgs <- renderPlot({
        plot.new()
      })
      validate(
        need(exists("spData"), "sparse data represention is not available. Please upload and load data first.")
      )
    }
    
    if(!exists("spwin")){
      output$message <- renderText({
        "spatial window is not available. Please upload and load data first."
      })
      output$imgs <- renderPlot({
        plot.new()
      })
      validate(
        need(exists("spwin"), "spatial window is not available. Please upload and load data first.")
      )
    }
    
    
    if(!is.numeric(input$mz)){
      output$message <- renderText({
        "The m/z value provided for MPM computation does not seem to be numeric!"
      })
      output$imgs <- renderPlot({
        plot.new()
      })
      validate(
        need(is.numeric(input$mz), "The m/z value provided for MPM computation does not seem to be numeric!")
      )
    }
    
    
    if(!spatstat.utils::check.in.range(input$mz, range(spData$mzAxis), FALSE)){
      output$message <- renderText({
        "The m/z value provided is outside of the mass range of the MSI data!"
      })
      output$imgs <- renderPlot({
        plot.new()
      })
      validate(
        need(spatstat.utils::check.in.range(input$mz, range(spData$mzAxis), FALSE), 
             "The m/z value provided is outside of the mass range of the MSI data!")
      )
    }
    
    
    
    
    
    # start loading
    withProgress(
      message="please wait",
      detail="Estimating molecular probabilitic map ...",
      value=0.2,{ 
        
        n <- 3
        
        # check if adjust is needed
        if(input$adjustmz){
          # use the nearest detectable m/z value
          querym <- spData$mzAxis[which.min(abs(spData$mzAxis - input$mz))]
          }
        else{
          querym <- input$mz
          }
        
        # search analyte
        incProgress(1/n, detail = paste("Searching for query mass .. "))
        sppMoi          <- moleculaR::searchAnalyte(m = querym, 
                                         fwhm = moleculaR::getFwhm(fwhmObj, querym), 
                                         spData = spData,
                                         spwin = spwin, 
                                         wMethod = "Gaussian")
        
        # another check
        if(sppMoi$n == 0){
          output$message <- renderText({
            "No signal has been detected for the provided m/z value! "
          })
          output$imgs <- renderPlot({
            plot.new()
          })
          validate(
            need(sppMoi$n > 0, 
                 "No signal has been detected for the provided m/z value! ")
          )
        }
        
        
        
        # compute MPM - default parameters
        incProgress(1/n, detail = paste("Estimating MPM .. "))
        probImg         <- moleculaR::probMap(sppMoi)
        
        # also compute an ion image for comparison
        sppIonImage      <- moleculaR::searchAnalyte(m = querym, 
                                          fwhm = moleculaR::getFwhm(fwhmObj, querym), 
                                          spData = spData, 
                                          spwin = spwin, 
                                          wMethod = "sum")
        
        
        # compute a raster image of the sppIonImage 
        ionImage        <- moleculaR::spp2im(sppIonImage)
        
        incProgress(1/n, detail = paste("Rendering plot .. "))
        
        output$message <- renderText({
          paste("The following m/z value was detected with the corresponding mass resolution: ", 
                paste("m/z ", sppMoi$metaData$mzVals, collapse = " "),
                paste("R", round(sppMoi$metaData$mzVals/moleculaR::getFwhm(fwhmObj, sppMoi$metaData$mzVals)), collapse = " "),
                sep = "\n")
        })
        
        output$imgs <- renderPlot({
          
          if(input$detailedPlotMPM){
            plot(probImg, what = "detailed",  ionImage = ionImage)
          } else {
            par(mfrow = c(1,2))
            plot(probImg, what = "MPM")
            moleculaR::plotImg(ionImage, main = "Ion Image")
          }
          
        })
        
        
      })
    
    
  })
  
  # observeEvent for CPPM
  observeEvent(input$renderCPPM,{
    
    
    # checks
    
    if(!exists("fwhmObj")){
      
      output$message <- renderText({
        "No FWHM data is available. Please upload and load data first."
      })
      
      output$imgs <- renderPlot({
        plot.new()
      })
      
      validate(
        need(exists("fwhmObj"), "No FWHM data is available. Please upload and load data first.")
      )
      
      }
 
    if(!exists("msData")){
      
      output$message <- renderText({
        "MSI data is not loaded. Please upload and load data first."
      })
      
      output$imgs <- renderPlot({
        plot.new()
      })
      
      validate(
        need(exists("msData"), "MSI data is not loaded. Please upload and load data first.")
      )
      
    }
    
   
    
    if(!exists("spData")){
      output$message <- renderText({
        "sparse data represention is not available. Please upload and load data first."
      })
      output$imgs <- renderPlot({
        plot.new()
      })
      validate(
        need(exists("spData"), "sparse data represention is not available. Please upload and load data first.")
      )
    }
    
    if(!exists("spwin")){
      output$message <- renderText({
        "spatial window is not available. Please upload and load data first."
      })
      output$imgs <- renderPlot({
        plot.new()
      })
      validate(
        need(exists("spwin"), "spatial window is not available. Please upload and load data first.")
      )
    }
    
    # convert the mz string into a numeric vector
    queryms <- as.numeric(strsplit(input$mzStr, ",")[[1]])
    
    
    if(any(is.na(queryms))){
      output$message <- renderText({
        paste0("There was a problem in the provided comma-separated m/z values!")
      })
      output$imgs <- renderPlot({
        plot.new()
      })
      validate(
        need(all(!is.na(queryms)), "There was a problem in the provided comma-separated m/z values!")
      )
    }
    
    
    # check if values are in range
    inRange <- sapply(queryms, FUN = function(i){
      spatstat.utils::check.in.range(i, range(spData$mzAxis), FALSE)
    })
    
    if(!all(inRange)){
      output$message <- renderText({
        paste0("One of the provided m/z values is outside of the mass range of the MSI data!")
      })
      output$imgs <- renderPlot({
        plot.new()
      })
      validate(
        need(all(inRange), "One of the provided m/z values is outside of the mass range of the MSI data!")
      )
    }
    
    
    
    
    # start loading
    withProgress(
      message="please wait",
      detail="Estimating collective projection probabilitic map ...",
      value=0.2,{ 
        
        n <- 5
        
        
        # search analyte
        incProgress(1/n, detail = paste("Searching for m/z values .. "))
        sppCPPM          <- moleculaR::batchAnalyteSearch(spData = spData, 
                                                         fwhmObj = fwhmObj, 
                                                         spwin = spwin, 
                                                         m = queryms)
        
        # another check
        if(sppCPPM$n == 0){
          output$message <- renderText({
            "No signal has been detected for all provided m/z value! "
          })
          output$imgs <- renderPlot({
            plot.new()
          })
          validate(
            need(sppCPPM$n > 0, "No signal has been detected for all provided m/z value! ")
          )
        }
        
        
        
        
        
        # standadize intensity
        incProgress(1/n, detail = paste("Standardizing intensities .. "))
        sppCPPM <- moleculaR::transformIntensity(sppCPPM, method = "z-score")
        
        # filter duplicates
        incProgress(1/n, detail = paste("Filtering duplicates .. "))
        sppCPPM <- moleculaR::filterDuplicates(sppCPPM)
        
        # compute MPM - default parameters
        incProgress(1/n, detail = paste("Estimating CPPM .. "))
        probImg         <- moleculaR::probMap(sppCPPM)
        
        
        incProgress(1/n, detail = paste("Rendering plot .. "))
        
        output$message <- renderText({
          paste("The following m/z values were detected with the corresponding mass resolution: ", 
                paste("m/z", sppCPPM$metaData$mzVals, collapse = " "),
                paste("R ", round(sppCPPM$metaData$mzVals/moleculaR::getFwhm(fwhmObj, sppCPPM$metaData$mzVals)), collapse = " "),
                sep = "\n")
        })
        
        output$imgs <- renderPlot({
          
          if(input$detailedPlotCPPM){
            plot(probImg, what = "detailed")
          } else {
            plot(probImg, what = "MPM")
          }
          
        })
        
        
      })
    
    
  })
  

  # help text
  # output$helpTxt <- renderUI({
  #   HTML(paste("","","",helpTxt1, helpTxt2, sep="<br/>"))
  # })
  
  # exit app gracefully
  observeEvent(input$exit, {
    if(input$exit > 0){
      stopApp("App Exited.")
    }
  })
  
  
  # observeEvent for loading exampple
  observeEvent(input$loadExample,{
    
    
    
    # start loading
    withProgress(
      message="please wait",
      detail="Loading and Processing Data...",
      value=0.2,{ 
        
        n <- 6
        
        incProgress(1/n, detail = paste("Load example data .. "))
        data("processed-example-Data", package = "moleculaR")
        
        
        # create sparse matrix representation
        incProgress(1/n, detail = paste("Creating sparse matrix representation .. "))
        spData <<- moleculaR::createSparseMat(x = msData)
        
        # create a spatial window
        incProgress(1/n, detail = paste("Creating a spatial window .. "))
        
        output$imgsExample <- renderPlot({
          par(mfrow = c(2, 1))
          spwin     <<- moleculaR::createSpatialWindow(pixelCoords = MALDIquant::coordinates(msData), 
                                                       clean = TRUE,
                                                       plot = TRUE)
          
          plot(fwhmObj)
        })
        
        
      })
    
    
  })
  
  # observeEvent for MPM
  observeEvent(input$renderMPMExample,{
    
    
    
    
    if(!is.numeric(input$mzExample)){
      output$messageExapmle <- renderText({
        "The m/z value provided for MPM computation does not seem to be numeric!"
      })
      output$imgsExample <- renderPlot({
        plot.new()
      })
      validate(
        need(is.numeric(input$mzExample), "The m/z value provided for MPM computation does not seem to be numeric!")
      )
    }
    
    
    if(!spatstat.utils::check.in.range(input$mzExample, range(spData$mzAxis), FALSE)){
      output$messageExample <- renderText({
        "The m/z value provided is outside of the mass range of the MSI data!"
      })
      output$imgsExample <- renderPlot({
        plot.new()
      })
      validate(
        need(spatstat.utils::check.in.range(input$mzExample, range(spData$mzAxis), FALSE), 
             "The m/z value provided is outside of the mass range of the MSI data!")
      )
    }
    
    
    
    
    
    # start loading
    withProgress(
      message="please wait",
      detail="Estimating molecular probabilitic map ...",
      value=0.2,{ 
        
        n <- 3
        
        # check if adjust is needed
        if(input$adjustmzExample){
          # use the nearest detectable m/z value
          querym <- spData$mzAxis[which.min(abs(spData$mzAxis - input$mzExample))]
        }
        else{
          querym <- input$mzExample
        }
        
        # search analyte
        incProgress(1/n, detail = paste("Searching for query mass .. "))
        sppMoi          <- moleculaR::searchAnalyte(m = querym, 
                                                    fwhm = moleculaR::getFwhm(fwhmObj, querym), 
                                                    spData = spData,
                                                    spwin = spwin, 
                                                    wMethod = "Gaussian")
        
        # another check
        if(sppMoi$n == 0){
          output$messageExample <- renderText({
            "No signal has been detected for the provided m/z value! "
          })
          output$imgsExample <- renderPlot({
            plot.new()
          })
          validate(
            need(sppMoi$n > 0, 
                 "No signal has been detected for the provided m/z value! ")
          )
        }
        
        
        
        # compute MPM - default parameters
        incProgress(1/n, detail = paste("Estimating MPM .. "))
        probImg         <- moleculaR::probMap(sppMoi)
        
        # also compute an ion image for comparison
        sppIonImage      <- moleculaR::searchAnalyte(m = querym, 
                                                     fwhm = moleculaR::getFwhm(fwhmObj, querym), 
                                                     spData = spData, 
                                                     spwin = spwin, 
                                                     wMethod = "sum")
        
        
        # compute a raster image of the sppIonImage 
        ionImage        <- moleculaR::spp2im(sppIonImage)
        
        incProgress(1/n, detail = paste("Rendering plot .. "))
        
        output$messageExample <- renderText({
          paste("The following m/z value was detected with the corresponding mass resolution: ", 
                paste("m/z ", sppMoi$metaData$mzVals, collapse = " "),
                paste("R", round(sppMoi$metaData$mzVals/moleculaR::getFwhm(fwhmObj, sppMoi$metaData$mzVals)), collapse = " "),
                sep = "\n")
        })
        
        output$imgsExample <- renderPlot({
          
          if(input$detailedPlotMPMExample){
            plot(probImg, what = "detailed",  ionImage = ionImage)
          } else {
            par(mfrow = c(1,2))
            plot(probImg, what = "MPM")
            moleculaR::plotImg(ionImage, main = "Ion Image")
          }
          
        })
        
        
      })
    
    
  })
  
  # observeEvent for CPPM
  observeEvent(input$renderCPPMExample,{
    
    
    # convert the mz string into a numeric vector
    queryms <- as.numeric(strsplit(input$mzStrExample, ",")[[1]])
    
    
    if(any(is.na(queryms))){
      output$messageExample <- renderText({
        paste0("There was a problem in the provided comma-separated m/z values!")
      })
      output$imgsExample <- renderPlot({
        plot.new()
      })
      validate(
        need(all(!is.na(queryms)), "There was a problem in the provided comma-separated m/z values!")
      )
    }
    
    
    # check if values are in range
    inRange <- sapply(queryms, FUN = function(i){
      spatstat.utils::check.in.range(i, range(spData$mzAxis), FALSE)
    })
    
    if(!all(inRange)){
      output$messageExample <- renderText({
        paste0("One of the provided m/z values is outside of the mass range of the MSI data!")
      })
      output$imgsExample <- renderPlot({
        plot.new()
      })
      validate(
        need(all(inRange), "One of the provided m/z values is outside of the mass range of the MSI data!")
      )
    }
    
    
    
    
    # start loading
    withProgress(
      message="please wait",
      detail="Estimating collective projection probabilitic map ...",
      value=0.2,{ 
        
        n <- 5
        
        
        # search analyte
        incProgress(1/n, detail = paste("Searching for m/z values .. "))
        sppCPPM          <- moleculaR::batchAnalyteSearch(spData = spData, 
                                                          fwhmObj = fwhmObj, 
                                                          spwin = spwin, 
                                                          m = queryms)
        
        # another check
        if(sppCPPM$n == 0){
          output$messageExample <- renderText({
            "No signal has been detected for all provided m/z value! "
          })
          output$imgsExample <- renderPlot({
            plot.new()
          })
          validate(
            need(sppCPPM$n > 0, "No signal has been detected for all provided m/z value! ")
          )
        }
        
        
        
        
        
        # standadize intensity
        incProgress(1/n, detail = paste("Standardizing intensities .. "))
        sppCPPM <- moleculaR::transformIntensity(sppCPPM, method = "z-score")
        
        # filter duplicates
        incProgress(1/n, detail = paste("Filtering duplicates .. "))
        sppCPPM <- moleculaR::filterDuplicates(sppCPPM)
        
        # compute MPM - default parameters
        incProgress(1/n, detail = paste("Estimating CPPM .. "))
        probImg         <- moleculaR::probMap(sppCPPM)
        
        
        incProgress(1/n, detail = paste("Rendering plot .. "))
        
        output$messageExample <- renderText({
          paste("The following m/z values were detected with the corresponding mass resolution: ", 
                paste("m/z", sppCPPM$metaData$mzVals, collapse = " "),
                paste("R ", round(sppCPPM$metaData$mzVals/moleculaR::getFwhm(fwhmObj, sppCPPM$metaData$mzVals)), collapse = " "),
                sep = "\n")
        })
        
        output$imgsExample <- renderPlot({
          
          if(input$detailedPlotCPPMExample){
            plot(probImg, what = "detailed")
          } else {
            plot(probImg, what = "MPM")
          }
          
        })
        
        
      })
    
    
  })
  
  # exit app gracefully
  observeEvent(input$exitExample, {
    if(input$exit > 0){
      stopApp("App Exited.")
    }
  })
  
  
}


# Run the application
shinyApp(ui = ui, server = server)
