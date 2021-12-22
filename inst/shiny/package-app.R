library(shiny)
library(shinyWidgets)
library(moleculaR)
library(shinythemes)
library(stringr)

###precomputing####

searchList           = initLipidSearch(swissdb = sldb)


#// initialize the swisslipids database

searchList_reactive <<- initLipidSearch(swissdb = sldb)


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
                                   column(12,
                                          HTML(paste0("<b>","Upload your Files","</b>")),
                                          fileInput("imzmlFile", NULL, multiple = TRUE, buttonLabel = "imzML & ibd File",
                                                    accept = c(".imzML",".ibd")
                                          )),
                                   column(12,
                                          fileInput("spectrFile", NULL, multiple = FALSE, buttonLabel = "Spectrum .tsv File",
                                                    accept = c(".tsv")
                                          )),
                                   column(12,
                                          fileInput("mtspcFile", NULL, multiple = FALSE, buttonLabel = "Metaspace .csv File",
                                                    accept = c(".csv")
                                          )),
                                ),
                                hr(),

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
                                          selectInput(inputId = "lipidSat", label = "Collective Projection Maps - Lipid Saturation", choices = c("sat", "mono-unsat.", "di-unsat", "poly-unsat")),
                                          actionButton(inputId = "go_lipid_sat", label = "Generate Plot",style='padding:6px; font-size:80%')

                                   ))

                                , width = 3),
                             mainPanel(
                                plotOutput("imgs", width = 600, height=900)

                             )

                          )),
                 tabPanel("About",
                          p("This is a web app that processes MSI data.
                            moleculaR is an open-source R package available at github.com/CeMOS/molecularR.",
                            style = "font-size:16px"))
)

####Backend ####

server <- function(input, output, session) {
   #increase maximum file uploade size to 48 GB
   options(shiny.maxRequestSize=48000*1024^2)
   # create empty reactive values
   rv <- reactiveValues(go_mz = list(),
                        go_lipid = list(searchList = searchList)
   )

   msData <- reactiveVal()
   spmat <- reactiveVal()
   msSpectr <- reactiveVal()
   fwhmObj <- reactiveVal()
   mtspc <- reactiveVal()
   msCoordinates <- reactiveVal()
   env <- environment()
   spwin <- reactiveVal()
   spwiny <- reactiveVal()
   .uniqueMass <- reactiveVal()
   isVerified <- reactiveVal()

   # reactive values for which image to output
   plot_output <- reactiveVal("initial")
   searchListcreated <- reactiveVal(FALSE)
   lysocreated <- reactiveVal(FALSE)

   observeEvent(input$adjustmz,{

      if(input$adjustmz==TRUE & plot_output()!="initial"){
         #use the nearest mz in metaspace
         rv$go_mz$mz <- as.numeric(mtspc$mz)[which.min(abs(as.numeric(mtspc$mz) - input$mz))]
      }
      else{
         rv$go_mz$mz <- input$mz
      }
   })

   observeEvent(input$mz,{

      mtspc_tmp = c()
      if(length(input$mtspcFile$datapath)!=0){
         mtspc_tmp <- read.csv(file = input$mtspcFile$datapath, skip = 2,header = TRUE, colClasses = "character")
      }

      if(input$adjustmz==TRUE & length(mtspc_tmp)!=0){
         #use the nearest mz in the dataset
         rv$go_mz$mz <- as.numeric(mtspc_tmp$mz)[which.min(abs(as.numeric(mtspc_tmp$mz) - input$mz))]
      }
      else{
         rv$go_mz$mz <- input$mz
      }



   })

   observeEvent(input$go_load, {
      withProgress(
         message="please wait",
         detail="Loading Data...",
         value=0.2,{
            n<-6

            if (length(input$imzmlFile$name)!=2){
               plot_output("loading_error")
               return()
            }

            if (sum(str_detect(input$imzmlFile$datapath, ".imzML"))!=1){
               plot_output("loading_error")
               return()
            }

            if (sum(str_detect(input$imzmlFile$datapath, ".ibd"))!=1){
               plot_output("loading_error")
               return()
            }

            if (str_detect(input$spectrFile$datapath, ".tsv")==FALSE){
               plot_output("loading_error")
               return()
            }

            if (str_detect(input$mtspcFile$datapath, ".csv")==FALSE){

               isVerified        <<- FALSE

            }

            if (str_detect(input$mtspcFile$datapath, ".csv")!=FALSE){

               isVerified        <<- TRUE

            }

            for(i in 1:length(input$imzmlFile$name)){
               file.copy(input$imzmlFile$datapath[i], paste0(tempdir(),"/", input$imzmlFile$name[i]))
            }
            file <- input$imzmlFile
            #find first imzml file in case more than one were uploaded
            imzmlidx <- which(str_detect(file$datapath, ".imzML"))[1]
            imzmlpath <- (paste0(tempdir(),"/", file$name[imzmlidx]))

            incProgress(1/n, detail = paste("Reading Data..."))

            if (length(imzmlidx)!=0){
               msData            <<- readCentrData(path = imzmlpath)
            }

            ### spectr file
            msSpectr          <<- readSingleSpect(isolate({input$spectrFile$datapath}))

            ### metaspacefile
            mtspc             <<- read.csv(file = input$mtspcFile$datapath, skip = 2,header = TRUE, colClasses = "character")

            if(exists("mtspc")){
               isVerified        <<- TRUE
            } else {
               isVerified        <<- FALSE
            }


            incProgress(1/n, detail = paste("Calculating fwhm..."))

            plot_output("show_fwhm")

            incProgress(1/n, detail = paste("Binning Data..."))

            msData            <<-  MALDIquant::binPeaks(msData,
                                                        tolerance = fwhmObj(400)/400, #focusing on lipids
                                                        method = "relaxed")

            msData            <<- filterPeaks(x = msData, minFreq = 0.01)

            incProgress(1/n, detail = paste("Creating sparse matrix..."))

            spwin <<- spatstat.geom::as.polygonal(spatstat.geom::owin(mask = as.data.frame(MALDIquant::coordinates(msData))))
            spData             <<- createSparseMat(x = msData)

            bwMethod          <<- "spAutoCor"

            gc()
         })
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
               #// Possible inputs: {

               adduct <- c("M-H", "M+H", "M+Na", "M+K") # <-- multiple choice possible, values =  c("M-H", "M+H", "M+Na", "M+K")
               confirmedOnly <- TRUE             # either TRUE of FALSE, whether to only include verified (confirmed) analytes.

               #}

               if(isVerified){
                  verifiedMasses = as.numeric(mtspc$mz)
               }else{
                  verifiedMasses = NA
               }


               searchList <<- batchLipidSearch(spData = spData, fwhmObj = fwhmObj, sldb = sldb,
                                              adduct = adduct, numCores = 4L,
                                              verifiedMasses = verifiedMasses,
                                              confirmedOnly = confirmedOnly, verbose = TRUE)

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
               searchList <<- batchLipidSearch(spData = spData, fwhmObj = fwhmObj, sldb = sldb,
                                               adduct = c("M+H", "M+Na", "M+K"), numCores = 4L, verifiedMasses = as.numeric(mtspc$mz),
                                               confirmedOnly = TRUE, verbose = TRUE)

               searchListcreated <<- reactiveVal(TRUE)
            }

            incProgress(1.5/n, detail = paste("Combining all lyso-GPLs into one SPP object"))

            if (lysocreated() != TRUE){

               ofInterest <- c("LPA(x:x)", "LPC(x:x)", "LPE(x:x)", "LPG(x:x)","LPI(x:x)", "LPS(x:x)",
                               "PA(x:x)", "PC(x:x)", "PE(x:x)","PG(x:x)", "PI(x:x)", "PS(x:x)")

               lysoGPLs <<- subsetAnalytes(searchList, lipidClass %in% ofInterest)

               lysocreated <<- reactiveVal(TRUE)

            }



            detectedAdducts <<- detectedSaturation <<- c("M+K", "M+Na", "M+H")

            sppList <<- setNames(vector("list", length(detectedAdducts)), detectedAdducts)

            # subsetting
            sppList[["M+K"]]     <<- subsetAnalytes(lysoGPLs, adduct == "M+K")
            sppList[["M+Na"]]    <<- subsetAnalytes(lysoGPLs, adduct == "M+Na")
            sppList[["M+H"]]     <<- subsetAnalytes(lysoGPLs, adduct == "M+H")

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
               searchList <<- batchLipidSearch(spData = spData, fwhmObj = fwhmObj, sldb = sldb,
                                               adduct = c("M+H", "M+Na", "M+K"), numCores = 4L, verifiedMasses = as.numeric(mtspc$mz),
                                               confirmedOnly = TRUE, verbose = TRUE)

               searchListcreated <<- reactiveVal(TRUE)
            }

            incProgress(1.5/n, detail = paste("Combining all lyso-GPLs into one SPP object"))

            if (lysocreated() != TRUE){

               ofInterest <- c("LPA(x:x)", "LPC(x:x)", "LPE(x:x)", "LPG(x:x)","LPI(x:x)", "LPS(x:x)",
                               "PA(x:x)", "PC(x:x)", "PE(x:x)","PG(x:x)", "PI(x:x)", "PS(x:x)")

               lysoGPLs <<- subsetAnalytes(searchList, lipidClass %in% ofInterest)

               lysocreated <<- reactiveVal(TRUE)

            }



            detectedSaturation <<- c("sat", "mono-unsat", "di-unsat", "poly-unsat")

            sppList <<- setNames(vector("list", length(detectedSaturation)), detectedSaturation)

            # subsetting
            sppList[["sat"]]   <<- subsetAnalytes(lysoGPLs, numDoubleBonds == 0)
            sppList[["mono-unsat"]]   <<- subsetAnalytes(lysoGPLs, numDoubleBonds == 1)
            sppList[["di-unsat"]]   <<- subsetAnalytes(lysoGPLs, numDoubleBonds == 2)
            sppList[["poly-unsat"]]   <<- subsetAnalytes(lysoGPLs, numDoubleBonds > 2)

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

            } else{ # if there are hits then proceed with MPM computations

               # compute rastered image of the sppIonImage
               ionImage        <- spatstat.geom::pixellate(rv$go_mz$hitsIonImage,
                                                      weights = rv$go_mz$hitsIonImage$marks$intensity,
                                                      W = spatstat.geom::as.mask(rv$go_mz$hitsIonImage$window,
                                                                            dimyx=c(diff(rv$go_mz$hitsIonImage$window$yrange),
                                                                                    diff(rv$go_mz$hitsIonImage$window$xrange))),
                                                      padzero = FALSE)


               # compute sppMoi (spatial point pattern of the analyte)
               sppMoi          <- searchAnalyte(m = rv$go_mz$mz_updated,
                                                fwhm = getFwhm(fwhmObj, rv$go_mz$mz_updated),
                                                spData = spData,
                                                wMethod = "Gaussian")




               #// compute MPM - default parameters
               probImg         <- probMap(sppMoi)

               txt  <- paste0("m/z ", round(rv$go_mz$mz_updated, 4), " ± ", round(getFwhm(fwhmObj, rv$go_mz$mz_updated), 4))
               plot(probImg, what = "detailed", analyte = txt, ionImage = ionImage)

               rm(probImg, txt, sppMoi, ionImage)

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

               probImg <- probMap(paHits, bwMethod = "scott", sqrtTansform = TRUE) # fixed arguments

               plot(probImg, what = "detailed", analyte = paste0(lipidClass_iso, " - n=", length(probImg$sppMoi$metaData$mzVals)))

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

            if(identical(sppList[[igroup]], NULL)) {



               par(mfrow = c(1, 1))

               #// empty window
               spatstat.geom::plot.owin(lysoGplsSumSpp$window,
                                   main = paste0("No insances of ", igroup, " were detected"),
                                   ylim = rev(lysoGplsSumSpp$window$yrange),
                                   box = FALSE)



            } else {

               probImg    = probMap(sppList[[igroup]], bwMethod = "scott", sqrtTansform = TRUE)

               if(probImg$sppMoi$n > 50000) {
                  cat("plotting ", format(probImg$sppMoi$n, big.mark = ","), " points - this takes time! \n")
               }
               plot(probImg, what = "detailed", analyte = paste0(igroup, " of ", lipidGroup, " - n=", length(probImg$sppMoi$metaData$mzVals)))

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

            if(identical(sppList[[igroup]], NULL)) {



               par(mfrow = c(1, 1))

               #// empty window
               spatstat.geom::plot.owin(lysoGplsSumSpp$window,
                                   main = paste0("No insances of ", igroup, " lipids were detected"),
                                   ylim = rev(lysoGplsSumSpp$window$yrange),
                                   box = FALSE)


            } else {

               probImg    = probMap(sppList[[igroup]], bwMethod = "scott", sqrtTansform = TRUE)

               if(probImg$sppMoi$n > 50000) {
                  cat("plotting ", format(probImg$sppMoi$n, big.mark = ","), " points - this takes time! \n")
               }
               plot(probImg, what = "detailed", analyte = paste0(igroup, " of ", lipidGroup, " - n=", length(probImg$sppMoi$metaData$mzVals)))


            }


         })






   })

   # plot for imgs
   output$imgs <- renderPlot({
      #// plotting
      if(plot_output()=="loading_error"){
         plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
         text(x = 0.5, y = 0.5, paste("Loading Error\n check if your uploaded .imzml matches your .ibd file\n as well if the\n .tsv and .csv are uploaded correct"),
              cex = 1.6, col = "black")
      }

      if(plot_output()=="initial"){
         plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
         text(x = 0.5, y = 0.5, paste("Initialize FWHM\n and then set your parameters\n to view and render the resulting plots.\n"),
              cex = 1.6, col = "black")
      }

      if(plot_output()=="show_fwhm"){

         fwhmObj           <<- estimateFwhm(s = msSpectr, plot = TRUE)

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
