library(shiny)
library(shinyWidgets)
library(moleculaR)
library(shinythemes)

data("processed-example-Data")
spData             = createSparseMat(x = msData)

#// initialize the swisslipids database
searchList           = initLipidSearch(swissdb = sldb)

####  Frontend ####
ui <- navbarPage(p("moleculaR: Spatial Probabilistic Mapping of Metabolites in Mass Spectrometry Imaging", HTML('&emsp;'), HTML('&emsp;')), theme = shinytheme("flatly"), selected="Main",
                 tabPanel("Main",
                          tags$head(
                             tags$style(HTML("hr {border-top: 1px solid #000000;}"))
                          ),
                          pageWithSidebar(
                             headerPanel(''),
                             sidebarPanel(
                                fluidRow(
                                   column(12,HTML(paste0("<b>","Peaks FWHM Estimation","</b>"))),
                                   column(12,
                                          actionButton(inputId = "go_load", label = "Initialize",style='padding:6px; font-size:80%'),
                                   )),
                                hr(),

                                fluidRow(
                                   column(12,HTML(paste0("<b>","Molecular Probability Maps","</b>"))),
                                   column(12,

                                          numericInput("mz", label = "", value = 788.5447),
                                          prettyCheckbox(
                                             inputId = "adjustmz", label = "Find Closest Detectable m/z", value = TRUE, icon = icon("check"), animation = "pulse"
                                          ),
                                          actionButton(inputId = "go_mz", label = "Generate Plot",style='padding:6px; font-size:80%'),

                                   )),
                                hr(),
                                fluidRow(
                                   column(12,
                                          selectInput(inputId = "lipidSpecies", label = "Collective Projection Maps - Lipid Classes", choices = searchList$allClasses, selected = "PI(x:x)"),
                                          actionButton(inputId = "go_lipid", label = "Generate Plot",style='padding:6px; font-size:80%')

                                   )),
                                hr(),
                                fluidRow(
                                   column(12,
                                          selectInput(inputId = "lipidIon", label = "Collective Projection Maps - Ion Milieu", choices = c("M+K", "M+Na", "M+H")),
                                          actionButton(inputId = "go_lipid_ion", label = "Generate Plot",style='padding:6px; font-size:80%')

                                   )),
                                hr(),
                                fluidRow(
                                   column(12,
                                          selectInput(inputId = "lipidSat", label = "Collective Projection Maps - Lipid Saturation", choices = c("sat", "mono-unsat", "di-unsat", "poly-unsat")),
                                          actionButton(inputId = "go_lipid_sat", label = "Generate Plot",style='padding:6px; font-size:80%')

                                   ))

                                , width = 3),
                             mainPanel(
                                plotOutput("imgs", width = 600, height=900)

                             )

                          )),
                 tabPanel("About",
                          p("This is an example web app with a preloaded sample MSI data of a wild-type Glioblastoma Multiform tissue sample.
                            MSI Data is restricted to Metaspace-confirmed lipids (SwissLipids DB) at 0.2 FDR in the positive ion mode.
                            moleculaR is an open-source R package available at github.com/CeMOS/molecularR.",
                            style = "font-size:16px"))
)


####Backend ####
server <- function(input, output, session) {

   # create empty reactive values
   rv <- reactiveValues(go_mz = list(),
                        go_lipid = list(searchList = searchList)
   )


   # reactive values for which image to output
   plot_output <- reactiveVal("initial")
   searchListcreated <- reactiveVal(FALSE)
   lysocreated <- reactiveVal(FALSE)

   observeEvent(input$adjustmz,{

      if(input$adjustmz==TRUE){
         #use the nearest mz in metaspace
         rv$go_mz$mz <- as.numeric(mtspc$mz)[which.min(abs(as.numeric(mtspc$mz) - input$mz))]
      }
      else{
         rv$go_mz$mz <- input$mz
      }
   })

   observeEvent(input$mz,{

      if(input$adjustmz==TRUE){
         #use the nearest mz in metaspace
         rv$go_mz$mz <- as.numeric(mtspc$mz)[which.min(abs(as.numeric(mtspc$mz) - input$mz))]
      }

      else{
         rv$go_mz$mz <- input$mz
      }



   })

   observeEvent(input$go_load, {
      plot_output("show_fwhm")
   })

   # routine for rv update in the m/z case
   observeEvent(input$go_mz, {

      withProgress(
         message="please wait",
         detail="Loading Data...",
         value=0.2,{
            n<-2

            updateNumericInput(session, "mz", value = round(rv$go_mz$mz,4))

            rv$go_mz$mz_updated = rv$go_mz$mz

            rv$go_mz$hitsIonImage              = searchAnalyte(m = rv$go_mz$mz_updated,
                                                               fwhm = getFwhm(fwhmObj, rv$go_mz$mz_updated),
                                                               spData = spData,
                                                               wMethod = "sum")

            })


      plot_output("mz")
   })

   # routine for calculations of lipid species
   observeEvent(input$go_lipid, {
      rv$go_lipid$lipidClass         <- input$lipidSpecies

      withProgress(
         message="please wait",
         detail="Batch lipid search is ongoing - this can take several minutes",
         value=0.1,{
            n <- 4
            # User needs to be notified that they have to wait

            if (searchListcreated() != TRUE){


               spwin     <- createSpatialWindow(pixelCoords = MALDIquant::coordinates(msData),
                                                   clean = TRUE,
                                                   plot = TRUE)

               searchList <<- batchLipidSearch(spData = spData, fwhmObj = fwhmObj, sldb = sldb, spwin = spwin,
                                               adduct = c("M+H", "M+Na", "M+K"), numCores = 1L, verifiedMasses = as.numeric(mtspc$mz),
                                               confirmedOnly = TRUE, verbose = TRUE)

               searchList <<- transformIntensity(searchList, method = "z-score")

               searchListcreated <<- reactiveVal(TRUE)
            }

            incProgress(1.5/n, detail = paste("Finished Lipid Hits"))

            plot_output("lipid")
         })
   })

   observeEvent(input$go_lipid_ion,{

      rv$go_lipid$lipidClass         <- input$lipidSpecies

      withProgress(
         message="please wait",
         detail="Batch lipid search is ongoing - this can take several minutes",
         value=0.1,{
            n <- 4
            # User needs to be notified that they have to wait

            if (searchListcreated() != TRUE){


                  spwin     <- createSpatialWindow(pixelCoords = MALDIquant::coordinates(msData),
                                                   clean = TRUE,
                                                   plot = TRUE)

                  searchList <<- batchLipidSearch(spData = spData, fwhmObj = fwhmObj, sldb = sldb, spwin = spwin,
                                                  adduct = c("M+H", "M+Na", "M+K"), numCores = 1L, verifiedMasses = as.numeric(mtspc$mz),
                                                  confirmedOnly = TRUE, verbose = TRUE)

                  searchList <<- transformIntensity(searchList, method = "z-score")

                  searchListcreated <<- reactiveVal(TRUE)

            }

            incProgress(1.5/n, detail = paste("Combining all lyso-GPLs into one SPP object"))

            if (lysocreated() != TRUE){

               ofInterest <- c("LPA(x:x)", "LPC(x:x)", "LPE(x:x)", "LPG(x:x)","LPI(x:x)", "LPS(x:x)",
                               "PA(x:x)", "PC(x:x)", "PE(x:x)","PG(x:x)", "PI(x:x)", "PS(x:x)")

               lysoGPLs <<- subsetAnalytes(searchList, lipidClass %in% ofInterest)

               lysocreated <<- reactiveVal(TRUE)
            }

            lipidGroup <<-"(lyso)GPLs"

            plot_output("lipid_ion")
         })

   })

   observeEvent(input$go_lipid_sat,{
      rv$go_lipid$lipidClass         <- input$lipidSpecies

      withProgress(
         message="please wait",
         detail="Batch lipid search is ongoing - this can take several minutes",
         value=0.1,{
            n <- 4
            # User needs to be notified that they have to wait

            if (searchListcreated() != TRUE){


                  spwin     <- createSpatialWindow(pixelCoords = MALDIquant::coordinates(msData),
                                                   clean = TRUE,
                                                   plot = TRUE)

                  searchList <<- batchLipidSearch(spData = spData, fwhmObj = fwhmObj, sldb = sldb, spwin = spwin,
                                                  adduct = c("M+H", "M+Na", "M+K"), numCores = 1L, verifiedMasses = as.numeric(mtspc$mz),
                                                  confirmedOnly = TRUE, verbose = TRUE)

                  searchList <<- transformIntensity(searchList, method = "z-score")

                  searchListcreated <<- reactiveVal(TRUE)
            }

            incProgress(1.5/n, detail = paste("Combining all lyso-GPLs into one SPP object"))

            if (lysocreated() != TRUE){
               ofInterest <- c("LPA(x:x)", "LPC(x:x)", "LPE(x:x)", "LPG(x:x)","LPI(x:x)", "LPS(x:x)",
                               "PA(x:x)", "PC(x:x)", "PE(x:x)","PG(x:x)", "PI(x:x)", "PS(x:x)")

               lysoGPLs <<- subsetAnalytes(searchList, lipidClass %in% ofInterest)

               lysocreated <<- reactiveVal(TRUE)
            }


            lipidGroup <<-"(lyso)GPLs"


            plot_output("lipid_sat")
         })


   })

   mz <- eventReactive(input$go_mz,{
      withProgress(
         message="please wait",
         detail="Calculating plot...",
         value=0.1,{
            n<-3

            incProgress(1/n, detail = paste("Generating plot...."))

            #// check if hits is empty
            if(rv$go_mz$hitsIonImage$n == 0)
            {

               par(mfrow = c(1, 1))
               #// image without masking
               spatstat.geom::plot.owin(rv$go_mz$hitsIonImage$window,
                                   main = paste0("No insances of m/z ", round(rv$go_mz$mz_updated, 4), " were detected"),
                                   ylim = rev(rv$go_mz$hitsIonImage$window$yrange),
                                   box = FALSE)

               # rm(rv$go_mz$hitsIonImage)

            } else{ # if there are hits then proceed with MPM computations

               # compute rastered image of the sppIonImage






               # compute sppMoi (spatial point pattern of the analyte)
               spwin     <- createSpatialWindow(pixelCoords = MALDIquant::coordinates(msData),
                                                   clean = TRUE,
                                                   plot = TRUE)

               sppMoi          <- searchAnalyte(m = rv$go_mz$mz_updated,
                                                fwhm = getFwhm(fwhmObj, rv$go_mz$mz_updated),
                                                spData = spData,
                                                spwin = spwin,
                                                wMethod = "Gaussian")




               #// compute MPM - default parameters
               probImg         <- probMap(sppMoi)

               plot(probImg, what = "detailed", ionImage = ionImage)

            }



         })
   })

   lipid <- eventReactive(input$go_lipid,{
      withProgress(
         message="please wait",
         detail="Generating first plot...",
         value=0.1,{
            n<-2

            lipidClass_iso <- isolate(input$lipidSpecies)

            # subset lipidHits
            paHits <- subsetAnalytes(searchList, lipidClass == lipidClass_iso)


            if(paHits$n==0) {

               par(mfrow = c(1, 1))
               #// empty window
               spwin = spatstat.geom::as.polygonal(spatstat.geom::owin(mask = as.data.frame(MALDIquant::coordinates(msData))))



               spatstat.geom::plot.owin(spwin,
                                   main = paste0("No insances of ", lipidClass_iso, " were detected"),
                                   ylim = rev(spwin$yrange),
                                   box = FALSE)


            } else {

               probImg <- probMap(paHits) # fixed arguments

               plot(probImg, what = "detailed")

               rm(probImg)


            }

         })

   })

   lipid_ion <- eventReactive(input$go_lipid_ion,{

      withProgress(
         message="please wait",
         detail="Plotting all lyso-GPLs points, this takes time",
         value=0.25,{

            n=5

            igroup              = isolate(input$lipidIon)

            spp_tmp <- subsetAnalytes(lysoGPLs, adduct == igroup)

            if(identical(spp_tmp, NULL)) {



               par(mfrow = c(1, 1))

               #// empty window
               spatstat.geom::plot.owin(lysoGPLs$window,
                                   main = paste0("No insances of ", igroup, " were detected"),
                                   ylim = rev(lysoGplsSumSpp$window$yrange),
                                   box = FALSE)



            } else {

               probImg    = probMap(spp_tmp)

               if(probImg$sppMoi$n > 50000) {
                  cat("plotting ", format(probImg$sppMoi$n, big.mark = ","), " points - this takes time! \n")
               }

               par(cex.lab = 2, cex.main = 2, cex.axis = 1.5)
               plot(probImg, what = "detailed")

               rm(probImg)

            }


            incProgress(1/n, detail = paste("Calculating Plots...."))



         })


   })

   lipid_sat <- eventReactive(input$go_lipid_sat,{

      withProgress(
         message="please wait",
         detail="Loading Data...",
         value=0.2,{

            n=4

            igroup = isolate(input$lipidSat)


            if(igroup=="sat"){
               spp_tmp <- subsetAnalytes(lysoGPLs, numDoubleBonds == 0)
            }            else if (igroup=="mono-unsat"){
               spp_tmp <- subsetAnalytes(lysoGPLs, numDoubleBonds == 1)
            }            else if (igroup=="di-unsat"){
               spp_tmp <- subsetAnalytes(lysoGPLs, numDoubleBonds == 2)
            }            else if(igroup=="poly-unsat"){
               spp_tmp <- subsetAnalytes(lysoGPLs, numDoubleBonds > 2)
            }



            if(identical(spp_tmp, NULL)) {



               par(mfrow = c(1, 1))

               #// empty window
               spatstat.geom::plot.owin(lysoGplsSumSpp$window,
                                   main = paste0("No insances of ", igroup, " lipids were detected"),
                                   ylim = rev(lysoGplsSumSpp$window$yrange),
                                   box = FALSE)


            } else {

               probImg    = probMap(spp_tmp)

               if(probImg$sppMoi$n > 50000) {
                  cat("plotting ", format(probImg$sppMoi$n, big.mark = ","), " points - this takes time! \n")
               }
               par(cex.lab = 2, cex.main = 2, cex.axis = 1.5)
               plot(probImg, what = "detailed")


            }


         })






   })

   # plot for imgs
   output$imgs <- renderPlot({
      #// plotting

      if(plot_output()=="initial"){
         plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
         text(x = 0.5, y = 0.5, paste("Initialize FWHM\n and then set your parameters\n to view and render the resulting plots.\n"),
              cex = 1.6, col = "black")
      }

      if(plot_output()=="show_fwhm"){

         plot(fwhmObj)

      }

      if(plot_output()=="mz"){

         mz()

      }

      if(plot_output()=="lipid"){

         lipid()

      }

      if(plot_output()=="lipid_ion"){

         lipid_ion()

      }

      if(plot_output()=="lipid_sat"){

         lipid_sat()

      }


   }, bg="transparent")


}

# Run the application
shinyApp(ui = ui, server = server)
