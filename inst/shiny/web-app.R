library(shiny)
library(shinyWidgets)
library(moleculaR)

spmat             = createSparseMat(x = msData)
fwhmFun           = fwhm$fwhmFun
searchList           = initLipidSearch(swissdb = sldb)

####  Frontend ####
ui <- fluidPage(

      titlePanel(title=htmlOutput("picture"),"MSI data Denis"),

      sidebarLayout(
            sidebarPanel(
                  fluidRow(
                        column(12,HTML(paste0("<b>","fwhm Correction","</b>"))),
                        column(12,
                               actionButton(inputId = "go_load", label = "Initialize")
                        ),
                        style = "border: 4px outset grey;"),
                  fluidRow(
                        column(12,
                               radioGroupButtons(
                                     inputId = "mz_bandwith",
                                     label = "m/z Bandwith Method",
                                     choices = c("scott","iterative")
                               ),

                               numericInput("mz", label = "Enter a m/z value", value = 496.33972),
                               prettyCheckbox(
                                     inputId = "adjustmz", label = "Find closeset m/z from Metaspace", value = TRUE, icon = icon("check"), animation = "pulse"
                               ),
                               actionButton(inputId = "go_mz", label = "m/z Plotting"),
                               style = "border: 4px outset grey;"
                        )),
                  fluidRow(
                        column(12,
                               radioGroupButtons(
                                     inputId = "lipid_bandwith",
                                     label = "Lipid Bandwith Method",
                                     choices = c("scott","iterative")
                               ),
                               selectInput(inputId = "lipidSpecies", label = "Lipid species", choices = searchList$allClasses, selected = "PI(x:x)"),
                               actionButton(inputId = "go_lipid", label = "Lipid Plotting"),
                               style = "border: 4px outset grey;"
                        ))

            ),

            mainPanel(
                  fluidRow(column(width = 12, plotOutput("imgs")),
                           column(width = 12, htmlOutput("conf_ids_msg"))
                  )
            )
      )

)


####Backend ####
server <- function(input, output, session) {

      # create empty reactive values
      rv <- reactiveValues(go_mz = list(),
                           go_lipid = list(searchList = searchList)
      )
      # reactive values for which image to output
      plot_output <- reactiveVal("initial")

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

      #image in header
      output$picture <-  renderText({
            c(
                  '<img src="',
                  "https://www.cemos.hs-mannheim.de/fileadmin/user_upload/dateien_global/sitelogo/cemos_logo_text_rgb_zw.svg",
                  '" width="100" height="40"',
                  '>'
            )
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

                        rv$go_mz$s                 = MALDIquant::msiSlices(x = msData, center = rv$go_mz$mz_updated ,
                                                                           tolerance = (fwhmFun(rv$go_mz$mz_updated ) / 2.355) * 3)

                        incProgress(1.5, detail = paste("Generated mz slices"))

                        #// output::molecular probability map
                        rv$go_mz$.uniqueMass       = as.numeric(spmat@Dimnames[[2]])


                        rv$go_mz$hits              = searchAnalyte(m = rv$go_mz$mz_updated ,
                                                                   fwhm = fwhmFun(rv$go_mz$mz_updated ),
                                                                   massAxis = rv$go_mz$.uniqueMass,
                                                                   spData = spmat,
                                                                   coords = MALDIquant::coordinates(msData))

                        #labelforplot
                        rv$go_mz$label <- rv$go_mz$mz_updated })


            plot_output("mz")
      })

      # routine for calculations of lipid species
      observeEvent(input$go_lipid, {
            rv$go_lipid$lipidClass         = input$lipidSpecies
            withProgress(
                  message="please wait",
                  detail="Loading Data...",
                  value=0.2,{
                        n<-4
                        .uniqueMass       = as.numeric(spmat@Dimnames[[2]])

                        incProgress(1.5/n, detail = paste("Loading Data..."))

                        lipidClass_iso <- isolate(rv$go_lipid$lipidClass)


                        searchList$lipidHits[[lipidClass_iso]] = parallel::mclapply(X = seq(1, nrow(searchList$swissList[[lipidClass_iso]])),
                                                                                            mc.cores = 5, FUN = function(i) {


                                                                                                  df            = data.frame(x = integer(0),
                                                                                                                             y = integer(0),
                                                                                                                             mass = numeric(0),
                                                                                                                             intensity = numeric(0),
                                                                                                                             adduct = character(0),
                                                                                                                             mode = character(0),
                                                                                                                             modeAdduct = character(0),
                                                                                                                             lipidID = character(0),
                                                                                                                             sumformula = character(0),
                                                                                                                             fullName = character(0),
                                                                                                                             abbrev = character(0),
                                                                                                                             numDoubleBonds = integer(0),
                                                                                                                             stringsAsFactors = F)

                                                                                                  msCoordinates = MALDIquant::coordinates(msData)

                                                                                                  #// protonated ----
                                                                                                  lipTmp        = searchList$swissList[[lipidClass_iso]]$`Exact m/z of [M+H]+`[i]

                                                                                                  if(!is.na(lipTmp)) { # for example there is no protonated version of this lipid species

                                                                                                        df     = rbind(df,
                                                                                                                       searchLipid(m = lipTmp,
                                                                                                                                   fwhm = fwhmFun(lipTmp),
                                                                                                                                   massAxis = .uniqueMass,
                                                                                                                                   spData = spmat,
                                                                                                                                   coords = msCoordinates,
                                                                                                                                   mtspc = mtspc, # <--
                                                                                                                                   adduct = "H1",  # <--
                                                                                                                                   mode = "positive",  # <--
                                                                                                                                   modeAdduct = "M+H",  # <--
                                                                                                                                   lipidID = searchList$swissList[[lipidClass_iso]]$`Lipid ID`[i],
                                                                                                                                   sumformula = searchList$swissList[[lipidClass_iso]]$`Formula (pH7.3)`[i],
                                                                                                                                   abbrev = searchList$swissList[[lipidClass_iso]]$`Abbreviation*`[i],
                                                                                                                                   numDoubleBonds = searchList$swissList[[lipidClass_iso]]$numDoubleBond[i]))


                                                                                                  }


                                                                                                  #// Na+ adduct ----
                                                                                                  lipTmp               = searchList$swissList[[lipidClass_iso]]$`Exact m/z of [M+Na]+`[i]

                                                                                                  if(!is.na(lipTmp)) { # for example there is no Na-adduct version


                                                                                                        df     = rbind(df,
                                                                                                                       searchLipid(m = lipTmp,
                                                                                                                                   fwhm = fwhmFun(lipTmp),
                                                                                                                                   massAxis = .uniqueMass,
                                                                                                                                   spData = spmat,
                                                                                                                                   coords = msCoordinates,
                                                                                                                                   mtspc = mtspc, # <--
                                                                                                                                   adduct = "Na1",  # <--
                                                                                                                                   mode = "positive",  # <--
                                                                                                                                   modeAdduct = "M+Na",  # <--
                                                                                                                                   lipidID = searchList$swissList[[lipidClass_iso]]$`Lipid ID`[i],
                                                                                                                                   sumformula = searchList$swissList[[lipidClass_iso]]$`Formula (pH7.3)`[i],
                                                                                                                                   abbrev = searchList$swissList[[lipidClass_iso]]$`Abbreviation*`[i],
                                                                                                                                   numDoubleBonds = searchList$swissList[[lipidClass_iso]]$numDoubleBond[i]))




                                                                                                  }


                                                                                                  #// K+ adduct ----
                                                                                                  lipTmp        = searchList$swissList[[lipidClass_iso]]$`Exact m/z of [M+K]+`[i]

                                                                                                  if(!is.na(lipTmp)) { # for example there is no Na-adduct version

                                                                                                        df     = rbind(df,
                                                                                                                       searchLipid(m = lipTmp,
                                                                                                                                   fwhm = fwhmFun(lipTmp),
                                                                                                                                   massAxis = .uniqueMass,
                                                                                                                                   spData = spmat,
                                                                                                                                   coords = msCoordinates,
                                                                                                                                   mtspc = mtspc, # <--
                                                                                                                                   adduct = "K1",  # <--
                                                                                                                                   mode = "positive",  # <--
                                                                                                                                   modeAdduct = "M+K",  # <--
                                                                                                                                   lipidID = searchList$swissList[[lipidClass_iso]]$`Lipid ID`[i],
                                                                                                                                   sumformula = searchList$swissList[[lipidClass_iso]]$`Formula (pH7.3)`[i],
                                                                                                                                   abbrev = searchList$swissList[[lipidClass_iso]]$`Abbreviation*`[i],
                                                                                                                                   numDoubleBonds = searchList$swissList[[lipidClass_iso]]$numDoubleBond[i]))




                                                                                                  }


                                                                                                  return(df)
                                                                                            })

                        incProgress(1/n, detail = paste("Updating Lipid Hits"))

                        searchList$lipidHits[[lipidClass_iso]] = do.call("rbind", searchList$lipidHits[[lipidClass_iso]])

                        incProgress(1.5/n, detail = paste("Finished Lipid Hits"))

                        #// added check
                        if(nrow(searchList$lipidHits[[lipidClass_iso]]) == 0) {

                              rv$go_lipid$detectionsExist = FALSE


                        } else {
                              rv$go_lipid$detectionsExist = TRUE
                        }


                        #// create the spatial point pattern
                        if(rv$go_lipid$detectionsExist){

                              rv$go_lipid$hitsppp              = spatstat::ppp(x = searchList$lipidHits[[lipidClass_iso]]$x,
                                                                               y = searchList$lipidHits[[lipidClass_iso]]$y,
                                                                               window = spwin,
                                                                               marks = searchList$lipidHits[[lipidClass_iso]][ , -c(1, 2)])

                        }

                        incProgress(1/n, detail = paste("Generating Output"))

                        if(rv$go_lipid$detectionsExist) {
                              rv$go_lipid$tmpm                 = unique(searchList$lipidHits[[lipidClass_iso]]$mass) # stores unique detected masses
                              rv$go_lipid$tmpc                 = unique(searchList$lipidHits[[lipidClass_iso]][c("mass", "confirmed", "modeAdduct")])
                        }

                        plot_output("lipid")
                  })
      })


      output$conf_ids_msg <- renderText({
            if(plot_output()=="lipid" ){
               if(rv$go_lipid$detectionsExist==TRUE){

                  paste("<b>",
                        "<br>",
                        "Confirmed(metaspace) | Detected(in dataset) | Total(swisslipids) = ",
                        length(which(rv$go_lipid$tmpc$confirmed)),
                        " | ",
                        length(rv$go_lipid$tmpm),
                        " | ",
                        nrow(rv$go_lipid$searchList$swissList[[rv$go_lipid$lipidClass]]) * 4,"<br>",
                        "<br>",
                        "<br>",
                        "M+H detected = ", length(which(rv$go_lipid$tmpc$modeAdduct == "M+H") ),
                        " | ",
                        "M+H confirmed = ", length(which(rv$go_lipid$tmpc$modeAdduct == "M+H" & rv$go_lipid$tmpc$confirmed)),
                        "<br>",
                        "M-H detected = ", length(which(rv$go_lipid$tmpc$modeAdduct == "M-H") ),
                        " | ",
                        "M-H confirmed = ", length(which(rv$go_lipid$tmpc$modeAdduct == "M-H" & rv$go_lipid$tmpc$confirmed)),
                        "<br>",
                        "M+Na detected = ", length(which(rv$go_lipid$tmpc$modeAdduct == "M+Na") ),
                        " | ",
                        "M+Na confirmed = ", length(which(rv$go_lipid$tmpc$modeAdduct == "M+Na" & rv$go_lipid$tmpc$confirmed == TRUE)),
                        "<br>",
                        "M+K detected = ", length(which(rv$go_lipid$tmpc$modeAdduct == "M+K") ),
                        " | ",
                        "M+K confirmed = ", length(which(rv$go_lipid$tmpc$modeAdduct == "M+K" & rv$go_lipid$tmpc$confirmed)),
                        "</br>")
               }
               else{
                  paste("<b>","Lipid Class: ", rv$go_lipid$lipidClass," not detected in given dataset", "</b>")
               }
            }
      })

      # plot for imgs
      output$imgs <- renderPlot({
            #// plotting

            if(plot_output()=="initial"){
                  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
                  text(x = 0.5, y = 0.5, paste("Initialize fwhm\n and then set your parameters\n to view the plots.\n"),
                       cex = 1.6, col = "black")
            }

            if(plot_output()=="show_fwhm"){

                  p           = fwhm$fwhmValues$peaks          #peaks of the single spectrum
                  r           = range(p)              # range of m/z
                  qp          = seq(r[1], r[2])       # query peaks

                  plot(x = p, y = fwhm$fwhmValues$fwhmValues,
                       main = "Estimated fwhm(m/z)", xlab = "m/z (Da)", ylab = "fwhm")
                  lines(x = qp, y = fwhmFun(qp), col = "green", lwd = 2)

            }

            if(plot_output()=="mz"){
                  par(mfrow=c(1,2))

                  withProgress(
                        message="please wait",
                        detail="Generating first plot...",
                        value=0.1,{
                              n<-2

                              if(sum(rv$go_mz$s, na.rm = T) == 0){
                                    #par(mfrow = c(1, 1))
                                    #// image without masking
                                    spatstat::plot.owin(spwin,
                                                        main = paste0("No instances of m/z ", round(rv$go_mz$mz_updated , 4), " were detected"),
                                                        ylim = rev(range(spwin$y)),
                                                        box = FALSE)
                              } else {

                                    MALDIquant::plotMsiSlice(x = rv$go_mz$s,
                                                             colRamp = colorRamp(viridis::viridis_pal(option = "inferno")(100)),
                                                             interpolate = F)

                              }

                              incProgress(1/n, detail = paste("Generating second plot...."))

                              #// check if hits is empty
                              if(nrow(rv$go_mz$hits) == 0)
                              {

                                    #par(mfrow = c(1, 1))
                                    #// image without masking
                                    spatstat::plot.owin(spwin,
                                                        main = paste0("No instances of m/z ", round(rv$go_mz$mz_updated , 4), " were detected"),
                                                        ylim = rev(range(spwin$y)),
                                                        box = FALSE)

                              } else {

                                    #par(mfrow = c(1, 1))


                                    #// create the spatial point pattern
                                    rv$go_mz$hitsppp              = spatstat::ppp(x = rv$go_mz$hits$x,
                                                                                  y = rv$go_mz$hits$y,
                                                                                  window = spwin,
                                                                                  marks = rv$go_mz$hits[ , -c(1, 2)])


                                    #// plotting

                                    probImg              = probMap(rv$go_mz$hitsppp, bwMethod = input$mz_bandwith)

                                    spatstat::plot.im(probImg$denspp,
                                                      main = paste0("m/z ", round(rv$go_mz$label,4), " +/- ", round((fwhmFun(rv$go_mz$label) / 2.355) * 3, 4)),
                                                      col = spatstat::colourmap(viridis::viridis_pal(option = "inferno")(100), range = range(probImg$denspp, na.rm = T)),
                                                      ylim = rev(range(spwin$y)),
                                                      box = FALSE)
                                    spatstat::plot.im(probImg$nonHotspotMask, col = rgb(0,0,0,0.3), add = TRUE)





                              }
                        })

            }

            if(plot_output()=="lipid"){

                  par(mfrow=c(1,2))

                  withProgress(
                        message="please wait",
                        detail="Generating first plot...",
                        value=0.1,{
                              n<-2

                              if(rv$go_lipid$detectionsExist) {



                                    probImg              = probMap(rv$go_lipid$hitsppp, bwMethod = input$lipid_bandwith)


                                    spatstat::plot.im(probImg$denspp,
                                                      main =  paste0(rv$go_lipid$lipidClass, " - all ", length(rv$go_lipid$tmpm), " detections"),
                                                      col = spatstat::colourmap(viridis::viridis_pal(option = "inferno")(100), range = range(probImg$denspp, na.rm = T)),
                                                      ylim = rev(range(spwin$y)),
                                                      box = FALSE)
                                    spatstat::plot.im(probImg$nonHotspotMask, col = rgb(0,0,0,0.3), add = TRUE)

                              }
                              else{
                                    par(mfrow = c(1, 2))

                                    #// empty window
                                    spatstat::plot.owin(spwin,
                                                        main = paste0("No instances of ", rv$go_lipid$lipidClass, " were detected"),
                                                        ylim = rev(range(spwin$y)),
                                                        box = FALSE)


                              }

                              incProgress(1/n, detail = paste("Generating second plot...."))

                              if(rv$go_lipid$detectionsExist) {

                                    subsetppp            = spatstat::subset.ppp(rv$go_lipid$hitsppp, confirmed == TRUE)

                                    #// check if empty
                                    if(subsetppp$n == 0) {

                                          #par(mfrow = c(1, 1))
                                          spatstat::plot.owin(spwin,
                                                              main = paste0(rv$go_lipid$lipidClass, " - ", "no confirmed detections"),
                                                              ylim = rev(range(spwin$y)),
                                                              box = FALSE)


                                    }

                                    else {

                                          probImg              = probMap(subsetppp,  bwMethod = input$lipid_bandwith)


                                          spatstat::plot.im(probImg$denspp,
                                                            main = paste0(rv$go_lipid$lipidClass, " - ", length(which(rv$go_lipid$tmpc$confirmed)), " confirmed detections"),
                                                            col = spatstat::colourmap(viridis::viridis_pal(option = "inferno")(100), range = range(probImg$denspp, na.rm = T)),
                                                            ylim = rev(range(spwin$y)),
                                                            box = FALSE)
                                          spatstat::plot.im(probImg$nonHotspotMask, col = rgb(0,0,0,0.3), add = TRUE)

                                    }
                              }
                        })







            }


      }, bg="transparent")


}

# Run the application
shinyApp(ui = ui, server = server)

