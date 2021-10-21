library(shiny)
library(shinyWidgets)
library(moleculaR)
library(shinythemes)

data("processed-example-Data")
spmat             = createSparseMat(x = msData)
fwhmFun           = fwhm$fwhmFun



#// load swiss lipid database -
# -- loaded -- #

#// filter sldb
Hidx = MALDIquant::match.closest(x = sldb$`Exact m/z of [M+H]+`,
                                 table = sort(as.numeric(mtspc$mz)),
                                 tolerance = fwhmFun(sldb$`Exact m/z of [M+H]+`))

Naidx = MALDIquant::match.closest(x = sldb$`Exact m/z of [M+Na]+`,
                                  table = sort(as.numeric(mtspc$mz)),
                                  tolerance = fwhmFun(sldb$`Exact m/z of [M+Na]+`))

kidx = MALDIquant::match.closest(x = sldb$`Exact m/z of [M+K]+`,
                                 table = sort(as.numeric(mtspc$mz)),
                                 tolerance = fwhmFun(sldb$`Exact m/z of [M+K]+`))



idxToKeep = c(which(!is.na(Hidx)), which(!is.na(Naidx)), which(!is.na(kidx)))

sldb = sldb[idxToKeep, ]

#// initialize the swisslipids database
searchList           = initLipidSearch(swissdb = sldb)


#// load metaspace identification list
# -- loaded --

#// check if these files exist
# -- no need do in this case --

#// find peaks
.uniqueMass       = as.numeric(spmat@Dimnames[[2]])



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
                                          selectInput(inputId = "lipidIon", label = "Collective Projection Maps - Ion Milieu", choices = c("all", "M+K", "M+Na", "M+H"), selected = "all"),
                                          actionButton(inputId = "go_lipid_ion", label = "Generate Plot",style='padding:6px; font-size:80%')

                                   )),
                                hr(),
                                fluidRow(
                                   column(12,
                                          selectInput(inputId = "lipidSat", label = "Collective Projection Maps - Lipid Saturation", choices = c("all", "sat.", "mono.", "di.", "poly."), selected = "all"),
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

   searchList_reactive <- reactiveVal()

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

            rv$go_mz$s                 = MALDIquant::msiSlices(x = msData, center = rv$go_mz$mz_updated ,
                                                               tolerance = (fwhmFun(rv$go_mz$mz_updated ) / 2.355) * 3)

            incProgress(1.5, detail = paste("Generated mz slices"))

            #// output::molecular probability map
            rv$go_mz$.uniqueMass       = as.numeric(spmat@Dimnames[[2]])


            rv$go_mz$hitsIonImage              = searchAnalyte(m = rv$go_mz$mz_updated ,
                                                               fwhm = fwhmFun(rv$go_mz$mz_updated ),
                                                               massAxis = rv$go_mz$.uniqueMass,
                                                               spData = spmat,
                                                               coords = MALDIquant::coordinates(msData),
                                                               wMethod = "sum")

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
         value=0.1,{
            n <- 4
            if (searchListcreated() != TRUE){
               n <- 30
               for(lipidClass in names(searchList$lipidHits)){
                  incProgress(1/n, detail = paste0("Loading Lipid: ", lipidClass))
                  #if(!(lipidClass %in% ofInterest)) {next}

                  searchList$lipidHits[[lipidClass]] = parallel::mclapply(X = seq(1, nrow(searchList$swissList[[lipidClass]])),
                                                                          mc.cores = 1, FUN = function(i) {


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
                                                                             lipTmp        = searchList$swissList[[lipidClass]]$`Exact m/z of [M+H]+`[i]

                                                                             if(!is.na(lipTmp)) { # for example there is no protonated version of this lipid species


                                                                                df     = rbind(df,
                                                                                               searchLipid(m = lipTmp,
                                                                                                           fwhm = fwhmFun(lipTmp),
                                                                                                           massAxis = .uniqueMass,
                                                                                                           spData = spmat,
                                                                                                           coords = msCoordinates,
                                                                                                           mtspc = mtspc, # <--
                                                                                                           confirmedOnly = TRUE,
                                                                                                           adduct = "H1",  # <--
                                                                                                           mode = "positive",  # <--
                                                                                                           modeAdduct = "M+H",  # <--
                                                                                                           lipidID = searchList$swissList[[lipidClass]]$`Lipid ID`[i],
                                                                                                           sumformula = searchList$swissList[[lipidClass]]$`Formula (pH7.3)`[i],
                                                                                                           abbrev = searchList$swissList[[lipidClass]]$`Abbreviation*`[i],
                                                                                                           numDoubleBonds = searchList$swissList[[lipidClass]]$numDoubleBond[i]))




                                                                             }


                                                                             #// Na+ adduct ----
                                                                             lipTmp               = searchList$swissList[[lipidClass]]$`Exact m/z of [M+Na]+`[i]

                                                                             if(!is.na(lipTmp)) { # for example there is no Na-adduct version



                                                                                df     = rbind(df,
                                                                                               searchLipid(m = lipTmp,
                                                                                                           fwhm = fwhmFun(lipTmp),
                                                                                                           massAxis = .uniqueMass,
                                                                                                           spData = spmat,
                                                                                                           coords = msCoordinates,
                                                                                                           mtspc = mtspc, # <--
                                                                                                           confirmedOnly = TRUE,
                                                                                                           adduct = "Na1",  # <--
                                                                                                           mode = "positive",  # <--
                                                                                                           modeAdduct = "M+Na",  # <--
                                                                                                           lipidID = searchList$swissList[[lipidClass]]$`Lipid ID`[i],
                                                                                                           sumformula = searchList$swissList[[lipidClass]]$`Formula (pH7.3)`[i],
                                                                                                           abbrev = searchList$swissList[[lipidClass]]$`Abbreviation*`[i],
                                                                                                           numDoubleBonds = searchList$swissList[[lipidClass]]$numDoubleBond[i]))




                                                                             }


                                                                             #// K+ adduct ----
                                                                             lipTmp        = searchList$swissList[[lipidClass]]$`Exact m/z of [M+K]+`[i]

                                                                             if(!is.na(lipTmp)) { # for example there is no Na-adduct version

                                                                                df     = rbind(df,
                                                                                               searchLipid(m = lipTmp,
                                                                                                           fwhm = fwhmFun(lipTmp),
                                                                                                           massAxis = .uniqueMass,
                                                                                                           spData = spmat,
                                                                                                           coords = msCoordinates,
                                                                                                           mtspc = mtspc, # <--
                                                                                                           confirmedOnly = TRUE,
                                                                                                           adduct = "K1",  # <--
                                                                                                           mode = "positive",  # <--
                                                                                                           modeAdduct = "M+K",  # <--
                                                                                                           lipidID = searchList$swissList[[lipidClass]]$`Lipid ID`[i],
                                                                                                           sumformula = searchList$swissList[[lipidClass]]$`Formula (pH7.3)`[i],
                                                                                                           abbrev = searchList$swissList[[lipidClass]]$`Abbreviation*`[i],
                                                                                                           numDoubleBonds = searchList$swissList[[lipidClass]]$numDoubleBond[i]))




                                                                             }


                                                                             return(df)
                                                                          })

                  #// merge
                  searchList$lipidHits[[lipidClass]] = do.call("rbind", searchList$lipidHits[[lipidClass]])
               }
               searchList_reactive <<- searchList
               searchListcreated <<- reactiveVal(TRUE)
            }

            incProgress(1.5/n, detail = paste("Finished Lipid Hits"))

            lipidClass_iso <- isolate(rv$go_lipid$lipidClass)

            rv$go_lipid$lipidclass_rows = nrow(searchList_reactive$lipidHits[[lipidClass_iso]])
            rv$go_lipid$lipidclass_x = searchList_reactive$lipidHits[[lipidClass_iso]]$x
            rv$go_lipid$lipidclass_y = searchList_reactive$lipidHits[[lipidClass_iso]]$y
            rv$go_lipid$lipidclass_marks = searchList_reactive$lipidHits[[lipidClass_iso]][ , -c(1, 2)]


            incProgress(1/n, detail = paste("Generating Output"))

            plot_output("lipid")
         })
   })

   observeEvent(input$go_lipid_ion,{

      plot_output("lipid_ion")
   })

   observeEvent(input$go_lipid_sat,{

      plot_output("lipid_sat")
   })

   mz <- eventReactive(input$go_mz,{
      withProgress(
         message="please wait",
         detail="Calculating plot...",
         value=0.1,{
            n<-3

            incProgress(1/n, detail = paste("Generating plot...."))

            #// check if hits is empty
            if(nrow(rv$go_mz$hitsIonImage) == 0)
            {

               par(mfrow = c(1, 1))
               #// image without masking
               spatstat::plot.owin(spwin,
                                   main = paste0("No instances of m/z ", round(rv$go_mz$mz_updated , 4), " were detected"),
                                   ylim = rev(range(spwin$y)),
                                   box = FALSE)


            } else {


               #// create the spatial point pattern
               rv$go_mz$sppMpm              = spatstat::ppp(x = rv$go_mz$hitsIonImage$x,
                                                            y = rv$go_mz$hitsIonImage$y,
                                                            window = spwin,
                                                            marks = rv$go_mz$hitsIonImage[ , -c(1, 2)])


               #// plotting

               rv$go_mz$sppMpm = spatstat::as.ppp(rv$go_mz$sppMpm)

               rv$go_mz$imgMpm  = spatstat::pixellate(rv$go_mz$sppMpm,
                                                      weights = rv$go_mz$sppMpm$marks$intensity,
                                                      W = spatstat::as.mask(spwin,dimyx=c(diff(spwin$yrange),diff(spwin$xrange))),
                                                      padzero = FALSE, savemap = FALSE)


               rv$go_mz$probImg              = probMap(rv$go_mz$sppMpm)

               par(mfrow = c(3, 2))

               # _________________________________________________ part one: CSR

               rv$mz$colfun        = spatstat::colourmap(col = spatstat::to.transparent((viridis::viridis_pal(option = "inferno")(100)), 0.7),
                                                         range = range(rv$go_mz$probImg$csrMoi$marks$intensity))

               spatstat::plot.ppp(rv$go_mz$probImg$csrMoi, use.marks = TRUE, which.marks = "intensity",
                                  ylim = rev(spwin$yrange),
                                  #cols = viridis::viridis_pal(option = "inferno")(100),
                                  #markscale = 0.000004,
                                  #zap = 0.0,
                                  #chars = 21,
                                  main = paste0("CSR at m/z ", round(rv$go_mz$mz_updated, 4), " ± ", round((fwhmFun(rv$go_mz$mz_updated) / 2.355) * 3, 4)),
                                  symap = spatstat::symbolmap(pch = 19,
                                                              cols = rv$mz$colfun,
                                                              size = 0.4,
                                                              range = range(rv$go_mz$probImg$csrMoi$marks$intensity))) # colors according to intensity

               # _________________________________________________ part two: SPP

               rv$mz$colfun        = spatstat::colourmap(col = spatstat::to.transparent((viridis::viridis_pal(option = "inferno")(100)), 0.7),
                                                         range = range(rv$go_mz$sppMpm$marks$intensity))

               spatstat::plot.ppp(rv$go_mz$sppMpm, use.marks = TRUE, which.marks = "intensity",
                                  ylim = rev(spwin$yrange),
                                  #cols = viridis::viridis_pal(option = "inferno")(100),
                                  #markscale = 0.000004,
                                  #zap = 0.0,
                                  #chars = 21,
                                  main = paste0("SPP at m/z ", round(rv$go_mz$mz_updated, 4), " ± ", round((fwhmFun(rv$go_mz$mz_updated) / 2.355) * 3, 4)),
                                  symap = spatstat::symbolmap(pch = 19,
                                                              cols = rv$mz$colfun,
                                                              size = 0.4,
                                                              range = range(rv$go_mz$sppMpm$marks$intensity))) # colors according to intensity



               # _________________________________________________ part three: rho of csr

               spatstat::plot.im(rv$go_mz$probImg$rhoCsr,
                                 main = expression(paste(rho["CSR"], "(x,y)")),
                                 col = spatstat::colourmap(viridis::viridis_pal(option = "inferno")(100), range = range(rv$go_mz$probImg$rhoCs, na.rm = T)),
                                 ylim = rev(range(spwin$y)),
                                 box = FALSE)


               # _________________________________________________ part four: rho of MOI

               spatstat::plot.im(rv$go_mz$probImg$rhoMoi,
                                 main = expression(paste(rho["MOI"], "(x,y)")),
                                 col = spatstat::colourmap(viridis::viridis_pal(option = "inferno")(100), range = range(rv$go_mz$probImg$rhoMoi, na.rm = T)),
                                 ylim = rev(range(spwin$y)),
                                 box = FALSE)



               # _________________________________________________ part five : regular ion image

               #// create the spatial point pattern
               rv$go_mz$sppIonImage            = spatstat::ppp(x = rv$go_mz$hitsIonImage$x,
                                                               y = rv$go_mz$hitsIonImage$y,
                                                               window = spwin,
                                                               marks = rv$go_mz$hitsIonImage[ , -c(1, 2)])

               rv$go_mz$sppIonImage            = spatstat::as.ppp(rv$go_mz$sppIonImage)


               # raster image - ion
               rv$go_mz$imgIon  = spatstat::pixellate(rv$go_mz$sppIonImage,
                                                      weights = rv$go_mz$sppIonImage$marks$intensity,
                                                      W = spatstat::as.mask(spwin,dimyx=c(diff(spwin$yrange),diff(spwin$xrange))),
                                                      padzero = FALSE, savemap = FALSE)

               spatstat::plot.im(rv$go_mz$imgIon,
                                 main = paste0("Ion Image at m/z ", round(rv$go_mz$mz_updated, 4), " ± ", round((fwhmFun(rv$go_mz$mz_updated) / 2.355) * 3, 4)),
                                 col = spatstat::colourmap(viridis::viridis_pal(option = "inferno")(100), range = range(rv$go_mz$imgIon, na.rm = T)),
                                 ylim = rev(range(spwin$y)),
                                 box = FALSE)


               # _________________________________________________ part six : MPM

               spatstat::plot.im(rv$go_mz$imgMpm,
                                 main = paste0("MPM at m/z ", round(rv$go_mz$mz_updated, 4)," ± ", round((fwhmFun(rv$go_mz$mz_updated) / 2.355) * 3, 4)),
                                 col = spatstat::colourmap(viridis::viridis_pal(option = "inferno")(100), range = range(rv$go_mz$imgMpm, na.rm = T)),
                                 ylim = rev(range(spwin$y)),
                                 box = FALSE)

               spatstat::plot.owin(rv$go_mz$probImg$hotspotpp$window, col = rgb(1,1,1,0.0), border = "white", lwd = 5,  add = TRUE)
               spatstat::plot.owin(rv$go_mz$probImg$hotspotpp$window, col = rgb(1,1,1,0.0), border = "red", lwd = 2.5, lty = "dashed",add = TRUE)

               spatstat::plot.owin(rv$go_mz$probImg$coldspotpp$window, col = rgb(1,1,1,0.0), border = "white", lwd = 5,  add = TRUE)
               spatstat::plot.owin(rv$go_mz$probImg$coldspotpp$window, col = rgb(1,1,1,0.0), border = "blue", lwd = 2.5, lty = "dashed",add = TRUE)

               legend("bottom", legend = c("Analyte Hotspot", "Analyte Coldspot"), lty = c("dashed"),
                      col = c("red", "blue"), bty = "n", horiz = TRUE)


            }
         })
   })

   lipid <- eventReactive(input$go_lipid,{
      withProgress(
         message="please wait",
         detail="Generating first plot...",
         value=0.1,{
            n<-2

            if(rv$go_lipid$lipidclass_rows==0) {



               par(mfrow = c(1, 1))

               #// empty window
               spatstat::plot.owin(spwin,
                                   main = paste0("No insances of ", rv$go_lipid$lipidClass, " were detected"),
                                   ylim = rev(range(spwin$y)),
                                   box = FALSE)

            }
            else{

               rv$go_lipid$sppCpm              = spatstat::ppp(x = rv$go_lipid$lipidclass_x,
                                                               y = rv$go_lipid$lipidclass_y,
                                                               window = spwin,
                                                               marks = rv$go_lipid$lipidclass_marks)


               rv$go_lipid$sppCpm              = spatstat::as.ppp(rv$go_lipid$sppCpm)

               rv$go$lipid$probImg              = probMap(rv$go_lipid$sppCpm, bwMethod = "scott", sqrtTansform = TRUE)

               par(mfrow = c(3, 2))

               # _________________________________________________ part one: CSR

               rv$go_lipid$colfun        = spatstat::colourmap(col = spatstat::to.transparent((viridis::viridis_pal(option = "inferno")(100)), 0.7),
                                                               range = range(rv$go$lipid$probImg$csrMoi$marks$intensity))

               spatstat::plot.ppp(rv$go$lipid$probImg$csrMoi, use.marks = TRUE, which.marks = "intensity",
                                  ylim = rev(spwin$yrange),
                                  #cols = viridis::viridis_pal(option = "inferno")(100),
                                  #markscale = 0.000004,
                                  #zap = 0.0,
                                  #chars = 21,
                                  main = paste0("CSR of ", rv$go_lipid$lipidClass),
                                  symap = spatstat::symbolmap(pch = 19,
                                                              cols = rv$go_lipid$colfun,
                                                              size = 0.4,
                                                              range = range(rv$go$lipid$probImg$csrMoi$marks$intensity))) # colors according to intensity

               # _________________________________________________ part two: SPP

               rv$go_lipid$colfun      = spatstat::colourmap(col = spatstat::to.transparent((viridis::viridis_pal(option = "inferno")(100)), 0.7),
                                                             range = range(rv$go$lipid$probImg$sppMoi$marks$intensity))

               spatstat::plot.ppp(rv$go$lipid$probImg$sppMoi, use.marks = TRUE, which.marks = "intensity",
                                  ylim = rev(spwin$yrange),
                                  #cols = viridis::viridis_pal(option = "inferno")(100),
                                  #markscale = 0.000004,
                                  #zap = 0.0,
                                  #chars = 21,
                                  main = paste0("SPP of ", rv$go_lipid$lipidClass),
                                  symap = spatstat::symbolmap(pch = 19,
                                                              cols = rv$go_lipid$colfun,
                                                              size = 0.4,
                                                              range = range(rv$go$lipid$probImg$sppMoi$marks$intensity))) # colors according to intensity



               # _________________________________________________ part three: rho of csr

               spatstat::plot.im(rv$go$lipid$probImg$rhoCsr,
                                 main = expression(paste(rho["CSR"], "(x,y)")),
                                 col = spatstat::colourmap(viridis::viridis_pal(option = "inferno")(100), range = range(rv$go$lipid$probImg$rhoCs, na.rm = T)),
                                 ylim = rev(range(spwin$y)),
                                 box = FALSE)


               # _________________________________________________ part four: rho of MOI

               spatstat::plot.im(rv$go$lipid$probImg$rhoMoi,
                                 main = expression(paste(rho["MOIs"], "(x,y)")),
                                 col = spatstat::colourmap(viridis::viridis_pal(option = "inferno")(100), range = range(rv$go$lipid$probImg$rhoMoi, na.rm = T)),
                                 ylim = rev(range(spwin$y)),
                                 box = FALSE)





               # _________________________________________________ part five : CPM

               spatstat::plot.im(rv$go$lipid$probImg$rhoMoi,
                                 main = paste0("CPM of ", rv$go_lipid$lipidClass),
                                 col = spatstat::colourmap(viridis::viridis_pal(option = "inferno")(100), range = range(rv$go$lipid$probImg$rhoMoi, na.rm = T)),
                                 ylim = rev(range(spwin$y)),
                                 box = FALSE)

               spatstat::plot.owin(rv$go$lipid$probImg$hotspotpp$window, col = rgb(1,1,1,0.0), border = "white", lwd = 5,  add = TRUE)
               spatstat::plot.owin(rv$go$lipid$probImg$hotspotpp$window, col = rgb(1,1,1,0.0), border = "red", lwd = 2.5, lty = "dashed",add = TRUE)

               spatstat::plot.owin(rv$go$lipid$probImg$coldspotpp$window, col = rgb(1,1,1,0.0), border = "white", lwd = 5,  add = TRUE)
               spatstat::plot.owin(rv$go$lipid$probImg$coldspotpp$window, col = rgb(1,1,1,0.0), border = "blue", lwd = 2.5, lty = "dashed",add = TRUE)

               legend("bottom", legend = c("Analyte Hotspot", "Analyte Coldspot"), lty = c("dashed"),
                      col = c("red", "blue"), bty = "n", horiz = TRUE)


            }
         })

   })

   lipid_ion <- eventReactive(input$go_lipid_ion,{

      withProgress(
         message="please wait",
         detail="Loading Data...",
         value=0.05,{

            n=5

            if (searchListcreated() != TRUE){
               n <- 30
               for(lipidClass in names(searchList$lipidHits)){
                  incProgress(1/n, detail = paste0("Loading Lipid: ", lipidClass))
                  #if(!(lipidClass %in% ofInterest)) {next}

                  searchList$lipidHits[[lipidClass]] = parallel::mclapply(X = seq(1, nrow(searchList$swissList[[lipidClass]])),
                                                                          mc.cores = 1, FUN = function(i) {


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
                                                                             lipTmp        = searchList$swissList[[lipidClass]]$`Exact m/z of [M+H]+`[i]

                                                                             if(!is.na(lipTmp)) { # for example there is no protonated version of this lipid species


                                                                                df     = rbind(df,
                                                                                               searchLipid(m = lipTmp,
                                                                                                           fwhm = fwhmFun(lipTmp),
                                                                                                           massAxis = .uniqueMass,
                                                                                                           spData = spmat,
                                                                                                           coords = msCoordinates,
                                                                                                           mtspc = mtspc, # <--
                                                                                                           confirmedOnly = TRUE,
                                                                                                           adduct = "H1",  # <--
                                                                                                           mode = "positive",  # <--
                                                                                                           modeAdduct = "M+H",  # <--
                                                                                                           lipidID = searchList$swissList[[lipidClass]]$`Lipid ID`[i],
                                                                                                           sumformula = searchList$swissList[[lipidClass]]$`Formula (pH7.3)`[i],
                                                                                                           abbrev = searchList$swissList[[lipidClass]]$`Abbreviation*`[i],
                                                                                                           numDoubleBonds = searchList$swissList[[lipidClass]]$numDoubleBond[i]))




                                                                             }


                                                                             #// Na+ adduct ----
                                                                             lipTmp               = searchList$swissList[[lipidClass]]$`Exact m/z of [M+Na]+`[i]

                                                                             if(!is.na(lipTmp)) { # for example there is no Na-adduct version



                                                                                df     = rbind(df,
                                                                                               searchLipid(m = lipTmp,
                                                                                                           fwhm = fwhmFun(lipTmp),
                                                                                                           massAxis = .uniqueMass,
                                                                                                           spData = spmat,
                                                                                                           coords = msCoordinates,
                                                                                                           mtspc = mtspc, # <--
                                                                                                           confirmedOnly = TRUE,
                                                                                                           adduct = "Na1",  # <--
                                                                                                           mode = "positive",  # <--
                                                                                                           modeAdduct = "M+Na",  # <--
                                                                                                           lipidID = searchList$swissList[[lipidClass]]$`Lipid ID`[i],
                                                                                                           sumformula = searchList$swissList[[lipidClass]]$`Formula (pH7.3)`[i],
                                                                                                           abbrev = searchList$swissList[[lipidClass]]$`Abbreviation*`[i],
                                                                                                           numDoubleBonds = searchList$swissList[[lipidClass]]$numDoubleBond[i]))




                                                                             }


                                                                             #// K+ adduct ----
                                                                             lipTmp        = searchList$swissList[[lipidClass]]$`Exact m/z of [M+K]+`[i]

                                                                             if(!is.na(lipTmp)) { # for example there is no Na-adduct version

                                                                                df     = rbind(df,
                                                                                               searchLipid(m = lipTmp,
                                                                                                           fwhm = fwhmFun(lipTmp),
                                                                                                           massAxis = .uniqueMass,
                                                                                                           spData = spmat,
                                                                                                           coords = msCoordinates,
                                                                                                           mtspc = mtspc, # <--
                                                                                                           confirmedOnly = TRUE,
                                                                                                           adduct = "K1",  # <--
                                                                                                           mode = "positive",  # <--
                                                                                                           modeAdduct = "M+K",  # <--
                                                                                                           lipidID = searchList$swissList[[lipidClass]]$`Lipid ID`[i],
                                                                                                           sumformula = searchList$swissList[[lipidClass]]$`Formula (pH7.3)`[i],
                                                                                                           abbrev = searchList$swissList[[lipidClass]]$`Abbreviation*`[i],
                                                                                                           numDoubleBonds = searchList$swissList[[lipidClass]]$numDoubleBond[i]))




                                                                             }


                                                                             return(df)
                                                                          })

                  #// merge
                  searchList$lipidHits[[lipidClass]] = do.call("rbind", searchList$lipidHits[[lipidClass]])
               }
               searchList_reactive <<- searchList
               searchListcreated <<- reactiveVal(TRUE)
            }


            incProgress(1/n, detail = paste("Calculating Plots...."))

            rv$go_lipid$ion$hitsList            = setNames(vector("list", 4), c("all", "M+K", "M+Na", "M+H"))

            rv$go_lipid$ion$hitsList[["all"]]   = do.call("rbind", searchList_reactive$lipidHits[c("PC(x:x)",
                                                                                                   "PE(x:x)",
                                                                                                   "PI(x:x)",
                                                                                                   "PS(x:x)",
                                                                                                   "LPC(x:x)",
                                                                                                   "LPE(x:x)",
                                                                                                   "LPI(x:x)",
                                                                                                   "LPS(x:x)")])

            #rv$go_lipid$ion$hitsList[["all"]]   = rv$go_lipid$ion$hitsList[["all"]][rv$go_lipid$ion$hitsList[["all"]]$confirmed == TRUE, ]

            rv$go_lipid$ion$hitsList[["all"]]   = spatstat::ppp(x = rv$go_lipid$ion$hitsList[["all"]]$x,
                                                                y = rv$go_lipid$ion$hitsList[["all"]]$y,
                                                                window = spwin,
                                                                marks = rv$go_lipid$ion$hitsList[["all"]][ , -c(1, 2)],
                                                                checkdup = FALSE)
            rv$go_lipid$ion$hitsList[["all"]]   = spatstat::as.ppp(rv$go_lipid$ion$hitsList[["all"]])

            rv$go_lipid$ion$hitsList[["M+K"]]     = spatstat::subset.ppp(rv$go_lipid$ion$hitsList[["all"]], modeAdduct == "M+K")

            rv$go_lipid$ion$hitsList[["M+Na"]]    = spatstat::subset.ppp(rv$go_lipid$ion$hitsList[["all"]], modeAdduct == "M+Na")
            rv$go_lipid$ion$hitsList[["M+H"]]     = spatstat::subset.ppp(rv$go_lipid$ion$hitsList[["all"]], modeAdduct == "M+H")

            lipidSpecies         = "alkali-(lyso)GPLs"


            igroup = input$lipidIon

            if(rv$go_lipid$ion$hitsList[[igroup]]$n == 0) {



               par(mfrow = c(1, 1))

               #// empty window
               spatstat::plot.owin(spwin,
                                   main = paste0("No insances of ", igroup, " were detected"),
                                   ylim = rev(range(spwin$y)),
                                   box = FALSE)


            } else {

               probImg    = probMap(rv$go_lipid$ion$hitsList[[igroup]], bwMethod = "scott", sqrtTansform = TRUE)




               par(mfrow = c(3, 2))

               incProgress(1/n, detail = paste("Generating Plots...."))

               # _________________________________________________ part one: CSR

               colfun        = spatstat::colourmap(col = spatstat::to.transparent((viridis::viridis_pal(option = "inferno")(100)), 0.7),
                                                   range = range(probImg$csrMoi$marks$intensity))

               spatstat::plot.ppp(probImg$csrMoi, use.marks = TRUE, which.marks = "intensity",
                                  ylim = rev(spwin$yrange),
                                  #cols = viridis::viridis_pal(option = "inferno")(100),
                                  #markscale = 0.000004,
                                  #zap = 0.0,
                                  #chars = 21,
                                  main = paste0("CSR of ", igroup),
                                  symap = spatstat::symbolmap(pch = 19,
                                                              cols = colfun,
                                                              size = 0.4,
                                                              range = range(probImg$csrMoi$marks$intensity))) # colors according to intensity

               # _________________________________________________ part two: SPP

               colfun        = spatstat::colourmap(col = spatstat::to.transparent((viridis::viridis_pal(option = "inferno")(100)), 0.7),
                                                   range = range(probImg$sppMoi$marks$intensity))

               spatstat::plot.ppp(probImg$sppMoi, use.marks = TRUE, which.marks = "intensity",
                                  ylim = rev(spwin$yrange),
                                  #cols = viridis::viridis_pal(option = "inferno")(100),
                                  #markscale = 0.000004,
                                  #zap = 0.0,
                                  #chars = 21,
                                  main = paste0("SPP of ", igroup),
                                  symap = spatstat::symbolmap(pch = 19,
                                                              cols = colfun,
                                                              size = 0.4,
                                                              range = range(probImg$sppMoi$marks$intensity))) # colors according to intensity



               # _________________________________________________ part three: rho of csr

               spatstat::plot.im(probImg$rhoCsr,
                                 main = expression(paste(rho["CSR"], "(x,y)")),
                                 col = spatstat::colourmap(viridis::viridis_pal(option = "inferno")(100), range = range(probImg$rhoCs, na.rm = T)),
                                 ylim = rev(range(spwin$y)),
                                 box = FALSE)


               # _________________________________________________ part four: rho of MOI

               spatstat::plot.im(probImg$rhoMoi,
                                 main = expression(paste(rho["MOIs"], "(x,y)")),
                                 col = spatstat::colourmap(viridis::viridis_pal(option = "inferno")(100), range = range(probImg$rhoMoi, na.rm = T)),
                                 ylim = rev(range(spwin$y)),
                                 box = FALSE)





               # _________________________________________________ part five : CPM

               spatstat::plot.im(probImg$rhoMoi,
                                 main = paste0("CPM of ", igroup),
                                 col = spatstat::colourmap(viridis::viridis_pal(option = "inferno")(100), range = range(probImg$rhoMoi, na.rm = T)),
                                 ylim = rev(range(spwin$y)),
                                 box = FALSE)

               spatstat::plot.owin(probImg$hotspotpp$window, col = rgb(1,1,1,0.0), border = "white", lwd = 5,  add = TRUE)
               spatstat::plot.owin(probImg$hotspotpp$window, col = rgb(1,1,1,0.0), border = "red", lwd = 2.5, lty = "dashed",add = TRUE)

               spatstat::plot.owin(probImg$coldspotpp$window, col = rgb(1,1,1,0.0), border = "white", lwd = 5,  add = TRUE)
               spatstat::plot.owin(probImg$coldspotpp$window, col = rgb(1,1,1,0.0), border = "blue", lwd = 2.5, lty = "dashed",add = TRUE)

               legend("bottom", legend = c("Analyte Hotspot", "Analyte Coldspot"), lty = c("dashed"),
                      col = c("red", "blue"), bty = "n", horiz = TRUE)

            }
         })


   })

   lipid_sat <- eventReactive(input$go_lipid_sat,{

      withProgress(
         message="please wait",
         detail="Loading Data...",
         value=0.2,{

            n=4

            if (searchListcreated() != TRUE){
               n <- 30
               for(lipidClass in names(searchList$lipidHits)){
                  incProgress(1/n, detail = paste0("Loading Lipid: ", lipidClass))
                  #if(!(lipidClass %in% ofInterest)) {next}

                  searchList$lipidHits[[lipidClass]] = parallel::mclapply(X = seq(1, nrow(searchList$swissList[[lipidClass]])),
                                                                          mc.cores = 1, FUN = function(i) {


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
                                                                             lipTmp        = searchList$swissList[[lipidClass]]$`Exact m/z of [M+H]+`[i]

                                                                             if(!is.na(lipTmp)) { # for example there is no protonated version of this lipid species


                                                                                df     = rbind(df,
                                                                                               searchLipid(m = lipTmp,
                                                                                                           fwhm = fwhmFun(lipTmp),
                                                                                                           massAxis = .uniqueMass,
                                                                                                           spData = spmat,
                                                                                                           coords = msCoordinates,
                                                                                                           mtspc = mtspc, # <--
                                                                                                           confirmedOnly = TRUE,
                                                                                                           adduct = "H1",  # <--
                                                                                                           mode = "positive",  # <--
                                                                                                           modeAdduct = "M+H",  # <--
                                                                                                           lipidID = searchList$swissList[[lipidClass]]$`Lipid ID`[i],
                                                                                                           sumformula = searchList$swissList[[lipidClass]]$`Formula (pH7.3)`[i],
                                                                                                           abbrev = searchList$swissList[[lipidClass]]$`Abbreviation*`[i],
                                                                                                           numDoubleBonds = searchList$swissList[[lipidClass]]$numDoubleBond[i]))




                                                                             }


                                                                             #// Na+ adduct ----
                                                                             lipTmp               = searchList$swissList[[lipidClass]]$`Exact m/z of [M+Na]+`[i]

                                                                             if(!is.na(lipTmp)) { # for example there is no Na-adduct version



                                                                                df     = rbind(df,
                                                                                               searchLipid(m = lipTmp,
                                                                                                           fwhm = fwhmFun(lipTmp),
                                                                                                           massAxis = .uniqueMass,
                                                                                                           spData = spmat,
                                                                                                           coords = msCoordinates,
                                                                                                           mtspc = mtspc, # <--
                                                                                                           confirmedOnly = TRUE,
                                                                                                           adduct = "Na1",  # <--
                                                                                                           mode = "positive",  # <--
                                                                                                           modeAdduct = "M+Na",  # <--
                                                                                                           lipidID = searchList$swissList[[lipidClass]]$`Lipid ID`[i],
                                                                                                           sumformula = searchList$swissList[[lipidClass]]$`Formula (pH7.3)`[i],
                                                                                                           abbrev = searchList$swissList[[lipidClass]]$`Abbreviation*`[i],
                                                                                                           numDoubleBonds = searchList$swissList[[lipidClass]]$numDoubleBond[i]))




                                                                             }


                                                                             #// K+ adduct ----
                                                                             lipTmp        = searchList$swissList[[lipidClass]]$`Exact m/z of [M+K]+`[i]

                                                                             if(!is.na(lipTmp)) { # for example there is no Na-adduct version

                                                                                df     = rbind(df,
                                                                                               searchLipid(m = lipTmp,
                                                                                                           fwhm = fwhmFun(lipTmp),
                                                                                                           massAxis = .uniqueMass,
                                                                                                           spData = spmat,
                                                                                                           coords = msCoordinates,
                                                                                                           mtspc = mtspc, # <--
                                                                                                           confirmedOnly = TRUE,
                                                                                                           adduct = "K1",  # <--
                                                                                                           mode = "positive",  # <--
                                                                                                           modeAdduct = "M+K",  # <--
                                                                                                           lipidID = searchList$swissList[[lipidClass]]$`Lipid ID`[i],
                                                                                                           sumformula = searchList$swissList[[lipidClass]]$`Formula (pH7.3)`[i],
                                                                                                           abbrev = searchList$swissList[[lipidClass]]$`Abbreviation*`[i],
                                                                                                           numDoubleBonds = searchList$swissList[[lipidClass]]$numDoubleBond[i]))




                                                                             }


                                                                             return(df)
                                                                          })

                  #// merge
                  searchList$lipidHits[[lipidClass]] = do.call("rbind", searchList$lipidHits[[lipidClass]])
               }
               searchList_reactive <<- searchList
               searchListcreated <<- reactiveVal(TRUE)
            }

            incProgress(1/n, detail = paste("Calculating Plots...."))

            rv$go_lipid$sat$hitsList[["all"]]   = do.call("rbind", searchList_reactive$lipidHits[c("PC(x:x)",
                                                                                                   "PE(x:x)",
                                                                                                   "PI(x:x)",
                                                                                                   "PS(x:x)",
                                                                                                   "LPC(x:x)",
                                                                                                   "LPE(x:x)",
                                                                                                   "LPI(x:x)",
                                                                                                   "LPS(x:x)")])


            rv$go_lipid$sat$hitsList[["all"]]   = spatstat::ppp(x = rv$go_lipid$sat$hitsList[["all"]]$x,
                                                                y = rv$go_lipid$sat$hitsList[["all"]]$y,
                                                                window = spwin,
                                                                marks = rv$go_lipid$sat$hitsList[["all"]][ , -c(1, 2)],
                                                                checkdup = FALSE)

            rv$go_lipid$sat$hitsList[["all"]]   = spatstat::as.ppp(rv$go_lipid$sat$hitsList[["all"]])

            rv$go_lipid$sat$hitsList[["sat."]]  = spatstat::subset.ppp(rv$go_lipid$sat$hitsList[["all"]], numDoubleBonds == 0L)
            rv$go_lipid$sat$hitsList[["mono."]]  = spatstat::subset.ppp(rv$go_lipid$sat$hitsList[["all"]], numDoubleBonds == 1L)
            rv$go_lipid$sat$hitsList[["di."]]  = spatstat::subset.ppp(rv$go_lipid$sat$hitsList[["all"]], numDoubleBonds == 2L)
            rv$go_lipid$sat$hitsList[["poly."]]  = spatstat::subset.ppp(rv$go_lipid$sat$hitsList[["all"]], numDoubleBonds > 2L)

            lipidSpecies         = "satrtn-(lyso)GPLs"


            igroup = input$lipidSat

            if(rv$go_lipid$sat$hitsList[[igroup]]$n == 0) {



               par(mfrow = c(1, 1))

               #// empty window
               spatstat::plot.owin(spwin,
                                   main = paste0("No insances of ", igroup, " were detected"),
                                   ylim = rev(range(spwin$y)),
                                   box = FALSE)


            } else {

               probImg    = probMap(rv$go_lipid$sat$hitsList[[igroup]], bwMethod = "scott", sqrtTansform = TRUE)



               par(mfrow = c(3, 2))

               incProgress(1/n, detail = paste("Generating Plots...."))

               # _________________________________________________ part one: CSR

               colfun        = spatstat::colourmap(col = spatstat::to.transparent((viridis::viridis_pal(option = "inferno")(100)), 0.7),
                                                   range = range(probImg$csrMoi$marks$intensity))

               spatstat::plot.ppp(probImg$csrMoi, use.marks = TRUE, which.marks = "intensity",
                                  ylim = rev(spwin$yrange),
                                  #cols = viridis::viridis_pal(option = "inferno")(100),
                                  #markscale = 0.000004,
                                  #zap = 0.0,
                                  #chars = 21,
                                  main = paste0("CSR of ", igroup),
                                  symap = spatstat::symbolmap(pch = 19,
                                                              cols = colfun,
                                                              size = 0.4,
                                                              range = range(probImg$csrMoi$marks$intensity))) # colors according to intensity

               # _________________________________________________ part two: SPP

               colfun        = spatstat::colourmap(col = spatstat::to.transparent((viridis::viridis_pal(option = "inferno")(100)), 0.7),
                                                   range = range(probImg$sppMoi$marks$intensity))

               spatstat::plot.ppp(probImg$sppMoi, use.marks = TRUE, which.marks = "intensity",
                                  ylim = rev(spwin$yrange),
                                  #cols = viridis::viridis_pal(option = "inferno")(100),
                                  #markscale = 0.000004,
                                  #zap = 0.0,
                                  #chars = 21,
                                  main = paste0("SPP of ", igroup),
                                  symap = spatstat::symbolmap(pch = 19,
                                                              cols = colfun,
                                                              size = 0.4,
                                                              range = range(probImg$sppMoi$marks$intensity))) # colors according to intensity



               # _________________________________________________ part three: rho of csr

               spatstat::plot.im(probImg$rhoCsr,
                                 main = expression(paste(rho["CSR"], "(x,y)")),
                                 col = spatstat::colourmap(viridis::viridis_pal(option = "inferno")(100), range = range(probImg$rhoCs, na.rm = T)),
                                 ylim = rev(range(spwin$y)),
                                 box = FALSE)


               # _________________________________________________ part four: rho of MOI

               spatstat::plot.im(probImg$rhoMoi,
                                 main = expression(paste(rho["MOIs"], "(x,y)")),
                                 col = spatstat::colourmap(viridis::viridis_pal(option = "inferno")(100), range = range(probImg$rhoMoi, na.rm = T)),
                                 ylim = rev(range(spwin$y)),
                                 box = FALSE)





               # _________________________________________________ part five : CPM

               spatstat::plot.im(probImg$rhoMoi,
                                 main = paste0("CPM of ", igroup),
                                 col = spatstat::colourmap(viridis::viridis_pal(option = "inferno")(100), range = range(probImg$rhoMoi, na.rm = T)),
                                 ylim = rev(range(spwin$y)),
                                 box = FALSE)

               spatstat::plot.owin(probImg$hotspotpp$window, col = rgb(1,1,1,0.0), border = "white", lwd = 5,  add = TRUE)
               spatstat::plot.owin(probImg$hotspotpp$window, col = rgb(1,1,1,0.0), border = "red", lwd = 2.5, lty = "dashed",add = TRUE)

               spatstat::plot.owin(probImg$coldspotpp$window, col = rgb(1,1,1,0.0), border = "white", lwd = 5,  add = TRUE)
               spatstat::plot.owin(probImg$coldspotpp$window, col = rgb(1,1,1,0.0), border = "blue", lwd = 2.5, lty = "dashed",add = TRUE)

               legend("bottom", legend = c("Analyte Hotspot", "Analyte Coldspot"), lty = c("dashed"),
                      col = c("red", "blue"), bty = "n", horiz = TRUE)
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

         p           = fwhm$fwhmValues$peaks          #peaks of the single spectrum
         r           = range(p)              # range of m/z
         qp          = seq(r[1], r[2])       # query peaks

         plot(x = p, y = fwhm$fwhmValues$fwhmValues,
              main = "Estimated fwhm(m/z)", xlab = "m/z (Da)", ylab = "fwhm")
         lines(x = qp, y = fwhmFun(qp), col = "green", lwd = 2)

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
