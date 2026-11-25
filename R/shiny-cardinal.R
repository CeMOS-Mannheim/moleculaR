# Shared module for both applications. Each instance owns its experiment,
# extraction cache and files; only portable numeric results leave the session.
.cardinalAnalysisUi <- function(id, source = c("upload", "example")) {
    source <- match.arg(source)
    ns <- shiny::NS(id)
    shiny::sidebarLayout(
        shiny::sidebarPanel(width = 4,
            if (source == "upload") shiny::fileInput(ns("imzmlFile"), "Centroided imzML and matching ibd", multiple = TRUE,
                accept = c(".imzML", ".imzml", ".ibd")) else
                shiny::helpText("Reduced centroided glioblastoma example: 32 x 32 coordinate crop. These data cannot estimate peak widths."),
            shiny::checkboxInput(ns("checkBinary"), "Strict imzML binary/metadata checks", FALSE),
            shiny::selectInput(ns("run"), "Experimental run", choices = character()),
            shiny::checkboxInput(ns("useRoi"), "Select a coordinate ROI", FALSE),
            shiny::conditionalPanel(sprintf("input['%s']", ns("useRoi")),
                shiny::numericInput(ns("xMin"), "Minimum x", 0), shiny::numericInput(ns("xMax"), "Maximum x", 100),
                shiny::numericInput(ns("yMin"), "Minimum y", 0), shiny::numericInput(ns("yMax"), "Maximum y", 100)),
            shiny::textInput(ns("pixelSpacing"), "Pixel spacing (optional x, y)", ""),
            shiny::radioButtons(ns("fwhmSource"), "FWHM source", c("Approximate instrument model" = "default", "Separate profile spectrum" = "spectrum")),
            shiny::conditionalPanel(sprintf("input['%s'] == 'default'", ns("fwhmSource")),
                shiny::selectInput(ns("fwhmModality"), "Instrument", c("FTICR", "TOF-reflector", "TIMS-TOF")),
                shiny::checkboxInput(ns("customMassRes"), "Custom resolving power", FALSE),
                shiny::conditionalPanel(sprintf("input['%s']", ns("customMassRes")),
                    shiny::numericInput(ns("fwhmMassRes"), "Resolving power", 123000, min = 1),
                    shiny::numericInput(ns("fwhmMass"), "Reference m/z", 249, min = 0.001))),
            shiny::conditionalPanel(sprintf("input['%s'] == 'spectrum'", ns("fwhmSource")),
                shiny::fileInput(ns("spectrFile"), "Separate profile spectrum TSV (m/z, intensity)", accept = ".tsv"),
                shiny::checkboxInput(ns("profileHeader"), "TSV has a header", FALSE),
                shiny::selectInput(ns("peakMethod"), "Cardinal peak picker", c("diff", "sd", "mad", "quantile", "filter", "cwt")),
                shiny::numericInput(ns("profileSnr"), "Profile SNR threshold", 5, min = 0)),
            shiny::selectInput(ns("normalization"), "Normalization", c("None" = "none", "TIC", "RMS", "median", "Internal standard" = "IS")),
            shiny::conditionalPanel(sprintf("input['%s'] != 'none'", ns("normalization")),
                shiny::selectInput(ns("normalizationScale"), "Normalization scale", c("Cardinal default" = "cardinal", "Unit" = "unit")),
                shiny::helpText("TIC uses only the features present in this file. Reduced data do not recover the original acquisition's full TIC.")),
            shiny::conditionalPanel(sprintf("input['%s'] == 'IS'", ns("normalization")),
                shiny::numericInput(ns("standardMz"), "Internal standard m/z", 400, min = 0.001),
                shiny::numericInput(ns("standardTol"), "Absolute m/z tolerance", 0.01, min = 0)),
            shiny::fileInput(ns("annotationFile"), "Verified masses CSV (optional mz column)", accept = ".csv"),
            shiny::helpText("Pre-filter annotations to your chosen FDR before upload. METASPACE CSVs with two metadata rows are also accepted."),
            shiny::actionButton(ns("loadData"), "Load and initialize", class = "btn-primary"),
            shiny::tags$hr(),
            shiny::selectInput(ns("wMethod"), "Intensity aggregation", c("Gaussian", "sum", "max", "mean")),
            shiny::checkboxInput(ns("confirmedOnly"), "Only annotation-confirmed results", source == "example"),
            shiny::checkboxInput(ns("standardize"), "Standardize each analyte in collective maps", TRUE),
            shiny::numericInput(ns("seed"), "Random seed", 1, min = 0, step = 1),
            shiny::selectInput(ns("bandwidthMode"), "KDE bandwidth", c("Automatic" = "auto", "Fixed" = "fixed")),
            shiny::conditionalPanel(sprintf("input['%s'] == 'fixed'", ns("bandwidthMode")),
                shiny::numericInput(ns("bandwidth"), "Bandwidth", 2, min = 0.001)),
            shiny::tabsetPanel(
                shiny::tabPanel("One analyte", shiny::numericInput(ns("mz"), "Theoretical m/z", 788.5447, min = 0.001),
                    shiny::actionButton(ns("renderMPM"), "Generate MPM")),
                shiny::tabPanel("Mass list", shiny::textAreaInput(ns("mzStr"), "Theoretical m/z values (comma or space separated)", "760.5851, 788.5447"),
                    shiny::actionButton(ns("renderCPPM"), "Generate CPPM")),
                shiny::tabPanel("Lipids", shiny::fileInput(ns("lipidFile"), "SwissLipids species TSV (optional)", accept = ".tsv"),
                    shiny::actionButton(ns("loadLipids"), "Load lipid database"),
                    shiny::selectInput(ns("lipidView"), "Group by", c("Lipid class" = "class", "(Lyso)GPL ion" = "ion", "(Lyso)GPL saturation" = "saturation")),
                    shiny::selectInput(ns("lipidClass"), "Lipid class", choices = character()),
                    shiny::selectInput(ns("lipidIon"), "(Lyso)GPL ion", c("[M+H]+", "[M+Na]+", "[M+K]+", "[M-H]-")),
                    shiny::selectInput(ns("lipidSat"), "Saturation", c("Saturated" = "sat", "Monounsaturated" = "mono", "Diunsaturated" = "di", "Polyunsaturated" = "poly")),
                    shiny::checkboxGroupInput(ns("adducts"), "Adducts", c("[M+H]+", "[M+Na]+", "[M+K]+", "[M-H]-"), selected = c("[M+H]+", "[M+Na]+", "[M+K]+")),
                    shiny::actionButton(ns("renderLipids"), "Generate lipid CPPM"))),
            shiny::tags$hr(), shiny::actionButton(ns("resetData"), "Clear session data")),
        shiny::mainPanel(width = 8,
            shiny::uiOutput(ns("status")),
            shiny::tags$details(shiny::tags$summary("Data notices"), shiny::verbatimTextOutput(ns("notices"))),
            shiny::selectInput(ns("plotView"), "Display", c("Probability map" = "MPM", "Aggregated intensity image" = "image", "Diagnostic views" = "detailed")),
            shiny::plotOutput(ns("imgs"), height = "780px"), shiny::textOutput(ns("summary")),
            shiny::tags$div(style = "display:none", shiny::textOutput(ns("hasResult"))),
            shiny::conditionalPanel(sprintf("output['%s'] == 'yes'", ns("hasResult")),
                shiny::downloadButton(ns("downloadRds"), "Analysis RDS"), shiny::downloadButton(ns("downloadCsv"), "Intensities CSV"),
                shiny::downloadButton(ns("downloadPng"), "Probability map PNG"))))
}

.shinyPlotAnalysis <- function(result, view = "MPM") {
    if (identical(view, "image")) {
        plotImg(result$ionImage, rescale = FALSE, main = paste(if (result$provenance$standardized) "Standardized" else "",
            result$provenance$aggregation, "aggregated intensity"))
    } else if (identical(view, "detailed")) {
        old <- graphics::par(mfrow = c(3, 2))
        on.exit(graphics::par(old))
        for (what in c("csrPointPattern", "analytePointPattern", "csrDensityImage", "analyteDensityImage", "kdeIntensitiesDistr", "MPM")) {
            graphics::plot(result$map, what = what)
        }
    } else graphics::plot(result$map, what = "MPM")
}

.cardinalAnalysisServer <- function(id, source = c("upload", "example"), examplePath = NULL) {
    source <- match.arg(source)
    if (is.null(examplePath)) examplePath <- system.file("extdata", "shiny-example", package = "moleculaR", mustWork = TRUE)
    shiny::moduleServer(id, function(input, output, session) {
        state <- shiny::reactiveValues(root = .shinySessionDirectory(), ready = FALSE, busy = FALSE,
            msData = NULL, fwhmObj = NULL, spwin = NULL, result = NULL, cache = .shinyCache(),
            sourcePath = NULL, sourceName = NULL, datasetDirectory = NULL, uploadStamp = NULL,
            normalizedPath = NULL, auxiliary = list(), runs = character(), coordinateBounds = NULL,
            appliedSettings = NULL, verifiedMasses = NA_real_, massLimits = NULL, lipids = NULL,
            lipidStamp = NULL, notices = character(), message = "Load and initialize data to begin.")
        clear <- function(removeFiles = FALSE) {
            state$ready <- FALSE
            state$msData <- state$result <- state$fwhmObj <- state$spwin <- NULL
            state$cache <- .shinyCache()
            state$appliedSettings <- state$coordinateBounds <- state$massLimits <- NULL
            state$verifiedMasses <- NA_real_
            state$lipids <- state$lipidStamp <- NULL
            if (removeFiles) {
                .shinyRemove(state$datasetDirectory, state$root)
                .shinyRemove(state$normalizedPath, state$root)
                for (aux in state$auxiliary) .shinyRemove(aux$path, state$root)
                state$sourcePath <- state$datasetDirectory <- state$normalizedPath <- state$uploadStamp <- NULL
                state$auxiliary <- list()
                state$runs <- character()
            }
            state$message <- "Session data cleared. Load and initialize data to begin."
        }
        session$onSessionEnded(function() shiny::isolate({
            clear()
            .shinyRemove(state$root, state$root, allowRoot = TRUE)
        }))
        task <- function(title, action, initializing = FALSE) {
            if (isTRUE(state$busy)) return(invisible(NULL))
            state$busy <- TRUE
            state$result <- NULL
            state$notices <- character()
            on.exit(state$busy <- FALSE)
            tryCatch(withCallingHandlers(shiny::withProgress(message = title, value = 0, {
                action(function(value, detail) shiny::setProgress(value = value, detail = detail))
            }), warning = function(w) {
                state$notices <- unique(c(state$notices, conditionMessage(w)))
                invokeRestart("muffleWarning")
            }), error = function(e) {
                state$result <- NULL
                if (initializing) state$ready <- FALSE
                state$message <- conditionMessage(e)
            })
            invisible(NULL)
        }
        updateClasses <- function() {
            classes <- sort(unique(state$lipids$lipidGroup))
            chosen <- .shinyValue(input$lipidClass, "")
            selected <- if (chosen %in% classes) chosen else if ("PC(x:x)" %in% classes) "PC(x:x)" else utils::head(classes, 1)
            shiny::updateSelectInput(session, "lipidClass", choices = classes, selected = selected)
        }
        shiny::observeEvent(input$loadData, {
            task("Initializing Cardinal data", function(progress) {
                .shinyInitialize(state, input, source, examplePath, progress)
                updateClasses()
            }, initializing = TRUE)
            shiny::updateSelectInput(session, "run", choices = state$runs,
                selected = if (length(state$runs) == 1L) state$runs else input$run)
            if (!is.null(state$coordinateBounds) && !isTRUE(input$useRoi)) {
                for (axis in c("x", "y")) {
                    shiny::updateNumericInput(session, paste0(axis, "Min"), value = state$coordinateBounds[[axis]][1])
                    shiny::updateNumericInput(session, paste0(axis, "Max"), value = state$coordinateBounds[[axis]][2])
                }
            }
        }, ignoreInit = TRUE)
        shiny::observe({
            applied <- state$appliedSettings
            current <- .shinyDataSettings(input, state$runs)
            upload <- .shinyUploadStamp(input$imzmlFile)
            if (!is.null(applied) && (!identical(current, applied) ||
                (source == "upload" && !identical(upload, state$uploadStamp)))) shiny::isolate({
                state$ready <- FALSE
                state$result <- NULL
                state$cache <- .shinyCache()
                state$message <- "Data or model settings changed. Load and initialize again."
            })
        })
        shiny::observeEvent(list(input$mz, input$mzStr, input$wMethod, input$confirmedOnly,
            input$standardize, input$seed, input$bandwidthMode, input$bandwidth,
            input$lipidClass, input$lipidView, input$lipidIon, input$lipidSat, input$adducts), {
            if (!is.null(state$result)) {
                state$result <- NULL
                state$message <- "Query or map settings changed. Generate the analysis again."
            }
        }, ignoreInit = TRUE)
        shiny::observeEvent(input$lipidFile, {
            state$lipids <- state$lipidStamp <- state$result <- NULL
        }, ignoreInit = TRUE, ignoreNULL = FALSE)
        shiny::observeEvent(input$loadLipids, {
            task("Loading lipid database", function(progress) {
                path <- .shinyAuxiliary(state, input$lipidFile, "lipids", "tsv")
                state$lipids <- if (is.null(path)) readRDS(file.path(examplePath,
                    if (source == "example") "lipids.rds" else "swisslipids-species.rds")) else
                    loadSwissDB(path)[, .shinyLipidColumns(), drop = FALSE]
                state$lipidStamp <- .shinyUploadStamp(input$lipidFile)
                updateClasses()
                state$message <- paste("Loaded", nrow(state$lipids), "lipid species.")
                progress(1, "Ready")
            })
        }, ignoreInit = TRUE)
        analyze <- function(kind) task("Calculating molecular probability map", function(progress) {
            state$result <- .shinyAnalyze(state, input, kind, progress)
            state$message <- "Analysis complete. Results and downloads reflect the current settings."
        })
        shiny::observeEvent(input$renderMPM, analyze("single"), ignoreInit = TRUE)
        shiny::observeEvent(input$renderCPPM, analyze("batch"), ignoreInit = TRUE)
        shiny::observeEvent(input$renderLipids, analyze("lipid"), ignoreInit = TRUE)
        shiny::observeEvent(input$resetData, clear(removeFiles = TRUE), ignoreInit = TRUE)
        output$status <- shiny::renderUI(shiny::tags$div(class = if (state$ready) "alert alert-success" else "alert alert-info", state$message))
        output$message <- shiny::renderText(state$message)
        output$notices <- shiny::renderText(paste(state$notices, collapse = "\n"))
        output$hasResult <- shiny::renderText(if (is.null(state$result)) "no" else "yes")
        shiny::outputOptions(output, "hasResult", suspendWhenHidden = FALSE)
        output$summary <- shiny::renderText({
            result <- state$result
            if (is.null(result)) return("")
            paste(result$provenance$aggregation, "aggregation;", result$provenance$acquiredPixels,
                "acquired pixels;", nrow(result$pointPattern$metaData), "analytes; standardization:", result$provenance$standardized)
        })
        output$imgs <- shiny::renderPlot({
            if (!is.null(state$result)) .shinyPlotAnalysis(state$result, .shinyValue(input$plotView, "MPM")) else if (state$ready) {
                old <- graphics::par(mfrow = c(1, 2))
                on.exit(graphics::par(old))
                graphics::plot(state$spwin, main = "Acquired tissue window")
                graphics::plot(state$fwhmObj)
            } else {
                graphics::plot.new()
                graphics::text(0.5, 0.5, "Load and initialize centroided MSI to begin.")
            }
        }, res = 100)
        output$downloadRds <- shiny::downloadHandler(function() "moleculaR-analysis.rds", function(file) .shinySaveAnalysis(state$result, file, "rds"))
        output$downloadCsv <- shiny::downloadHandler(function() "moleculaR-intensities.csv", function(file) .shinySaveAnalysis(state$result, file, "csv"))
        output$downloadPng <- shiny::downloadHandler(function() "moleculaR-probability-map.png", function(file) .shinySaveAnalysis(state$result, file, "png"))
        list(state = state, save = function(file, format = "rds") shiny::isolate(.shinySaveAnalysis(state$result, file, format)))
    })
}

.cardinalAppUi <- function(exampleOnly = FALSE) {
    if (is.null(getOption("shiny.maxRequestSize"))) options(shiny.maxRequestSize = 48 * 1024^3)
    tabs <- list()
    if (!exampleOnly) tabs <- c(tabs, list(shiny::tabPanel("Main", .cardinalAnalysisUi("main", "upload"))))
    tabs <- c(tabs, list(shiny::tabPanel("Example", .cardinalAnalysisUi("example", "example")),
        shiny::tabPanel("Help", shiny::tags$div(style = "max-width:850px;margin:30px;",
            shiny::tags$h3("Cardinal MSI analysis"),
            shiny::tags$p("Upload a matching centroided imzML/ibd pair, select one run and optional ROI, and initialize. Spectra stay file-backed; initialization never peak-picks MSI. The example is a reduced centroided crop."),
            shiny::tags$p("Choose an approximate instrument FWHM model or upload a separate profile spectrum TSV for measured widths. Never estimate widths from reduced centroided MSI."),
            shiny::tags$p("Enter theoretical masses without snapping them to an acquired peak. Choose aggregation, confirmation and map settings, then generate an MPM or collective map. Changing data or FWHM settings requires initialization again."),
            shiny::tags$p("The example includes verified annotations and a small lipid database. Main can use the packaged SwissLipids database or an uploaded species TSV. Duplicate lipid annotations remain separate analytes."),
            shiny::tags$p("Downloads contain portable maps, point intensities and provenance. RDS files do not contain the MSI experiment; the copied source and normalization files are removed when the session ends.")))))
    do.call(shiny::navbarPage, c(list(title = "moleculaR", header = shiny::tags$h2(style = "margin-left:20px;", "Molecular probability maps with Cardinal")), tabs))
}

.cardinalAppServer <- function(input, output, session, exampleOnly = FALSE) {
    modules <- list()
    if (!exampleOnly) modules$main <- .cardinalAnalysisServer("main", "upload")
    modules$example <- .cardinalAnalysisServer("example", "example")
    invisible(modules)
}
