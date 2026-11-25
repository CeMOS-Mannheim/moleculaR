# The development library has a newer-build startup notice; keep it separate
# from warnings produced by the scientific/session code under test.
invisible(suppressWarnings(suppressPackageStartupMessages(require("shiny", quietly = TRUE, character.only = TRUE))))

shiny_input <- function(upload = NULL, ...) {
    utils::modifyList(list(imzmlFile = upload, run = "", useRoi = FALSE, pixelSpacing = "1",
        checkBinary = FALSE, normalization = "none", normalizationScale = "cardinal",
        fwhmSource = "default", fwhmModality = "TOF-reflector", customMassRes = TRUE,
        fwhmMassRes = 1000, fwhmMass = 100, profileHeader = FALSE, peakMethod = "diff", profileSnr = 5,
        mz = 100.0005, mzStr = "100, 99.9, 100", wMethod = "Gaussian", confirmedOnly = FALSE,
        standardize = FALSE, seed = 1, bandwidthMode = "fixed", bandwidth = 1,
        lipidView = "class", lipidClass = "PC(x:x)", lipidIon = "[M-H]-", lipidSat = "mono", adducts = "[M-H]-"), list(...))
}

shiny_upload <- function(path) data.frame(name = basename(c(path, sub("[.]imzML$", ".ibd", path))),
    datapath = c(path, sub("[.]imzML$", ".ibd", path)))

shiny_state <- function() {
    state <- new.env(parent = emptyenv())
    state$root <- .shinySessionDirectory()
    state$auxiliary <- list()
    state$cache <- .shinyCache()
    state
}

shiny_fields <- function(pattern) {
    list(x = pattern$x, y = pattern$y, marks = pattern$marks$intensity,
        metadata = pattern$metaData[, setdiff(names(pattern$metaData), "idx"), drop = FALSE], window = pattern$window)
}

test_that("Shiny uploads, auxiliary files and cleanup are owned and bounded", {
    root <- .shinySessionDirectory()
    on.exit(.shinyRemove(root, root, TRUE))
    path <- write_msi_pair(make_msi_arrays(), tempfile(fileext = ".imzML"))
    original <- tools::md5sum(shiny_upload(path)$datapath)
    pair <- .shinyCopyPair(shiny_upload(path), root)
    expect_true(all(file.exists(c(pair$path, sub("imzML$", "ibd", pair$path)))))
    expect_identical(unname(tools::md5sum(shiny_upload(path)$datapath)), unname(original))
    expect_error(.shinyCopyPair(shiny_upload(path)[1, ], root), "matching")
    broken <- shiny_upload(path)
    broken$name[2] <- "other.ibd"
    expect_error(.shinyCopyPair(broken, root), "same stem")
    broken$name <- c("a.ibd", "a.ibd")
    expect_error(.shinyCopyPair(broken, root), "exactly one")
    expect_error(.shinyRemove(path, root), "outside")
    expect_error(.shinyRemove(tempdir(), root), "outside")
    expect_error(.shinyRemove(root, root), "outside")
    expect_true(file.exists(path))
    state <- shiny_state()
    on.exit(.shinyRemove(state$root, state$root, TRUE), add = TRUE)
    profile <- tempfile(fileext = ".tsv")
    writeLines("100\t2", profile)
    upload <- data.frame(name = "profile.tsv", datapath = profile)
    copied <- .shinyAuxiliary(state, upload, "profile", "tsv")
    unlink(profile)
    expect_identical(.shinyAuxiliary(state, upload, "profile", "tsv"), copied)
    expect_true(file.exists(copied))
    expect_error(.shinyAuxiliary(state, upload, "annotations", "csv"), "one csv")
    .shinyRemove(pair$directory, root)
    expect_false(file.exists(pair$path))
    expect_true(all(file.exists(shiny_upload(path)$datapath)))
})

test_that("Shiny input parsing and mass bounds do not construct a global intensity representation", {
    expect_equal(.shinyNumbers("100, 101;  102\n103"), 100:103)
    for (bad in c("", "100,,101", "100;", "abc", "NaN", "Inf", "-1", "0")) expect_error(.shinyNumbers(bad))
    annotations <- tempfile(fileext = ".csv")
    writeLines(c("mz,fdr", "100,0.05", "101,0.1"), annotations)
    expect_equal(.shinyAnnotations(annotations), c(100, 101))
    writeLines(c("# METASPACE metadata", "# Export metadata", "mz,fdr", "100,0.05"), annotations)
    expect_equal(.shinyAnnotations(annotations), 100)
    writeLines(c("mz", "-1"), annotations)
    expect_error(.shinyAnnotations(annotations), "positive")
    writeLines(c("mass", "100"), annotations)
    expect_error(.shinyAnnotations(annotations), "column named mz")
    testthat::local_mocked_bindings(.msiSpectrumVectors = function(...) stop("Read intensities"),
        createSparseMat = function(...) stop("Legacy matrix"), .package = "moleculaR")
    expect_equal(.shinyMassLimits(make_msi_arrays()), c(99.9, 101))
    expect_equal(.shinyMassLimits(make_msi_experiment()), c(99.9, 101))
    expect_error(.shinyMassLimits(make_msi_arrays(rep(list(numeric()), 3), rep(list(numeric()), 3))), "no valid")
})

test_that("Shiny initialization keeps file-backed data and rereads the original for normalization", {
    state <- shiny_state()
    on.exit(.shinyRemove(state$root, state$root, TRUE))
    path <- write_msi_pair(fixture_cardinal(), tempfile(fileext = ".imzML"))
    input <- shiny_input(shiny_upload(path))
    original <- tools::md5sum(shiny_upload(path)$datapath)
    testthat::local_mocked_bindings(readCentrData = function(...) stop("Legacy import"),
        createSparseMat = function(...) stop("Legacy matrix"), .package = "moleculaR")
    expect_true(.shinyInitialize(state, input, "upload", NULL))
    expect_true(state$ready)
    expect_s4_class(state$msData, "MSImagingArrays")
    expect_true(methods::is(Cardinal::intensity(state$msData), "matter"))
    expect_equal(.msiNumberOfPixels(state$msData), nrow(legacy_fixture()$coordinates))
    expect_equal(getFwhm(state$fwhmObj, 100), 0.1)
    expect_equal(.msiCoordinates(state$msData), legacy_fixture()$coordinates, ignore_attr = TRUE)
    input$normalization <- "TIC"
    input$normalizationScale <- "unit"
    .shinyInitialize(state, input, "upload", NULL)
    first <- .msiSpectrumVectors(state$msData, 1)$intensity
    derived <- state$normalizedPath
    expect_true(file.exists(derived))
    expect_equal(first, .msiSpectrumVectors(fixture_cardinal(), 1)$intensity /
        sum(.msiSpectrumVectors(fixture_cardinal(), 1)$intensity))
    .shinyInitialize(state, input, "upload", NULL)
    expect_equal(.msiSpectrumVectors(state$msData, 1)$intensity, first)
    expect_false(file.exists(derived))
    expect_identical(unname(tools::md5sum(shiny_upload(path)$datapath)), unname(original))
    unlink(shiny_upload(path)$datapath)
    input$normalization <- "none"
    .shinyInitialize(state, input, "upload", NULL)
    expect_equal(.msiSpectrumVectors(state$msData, 1)$intensity, .msiSpectrumVectors(fixture_cardinal(), 1)$intensity)
    expect_null(state$normalizedPath)
    testthat::local_mocked_bindings(normalizeIntensity = function(x, method, ..., outfile) {
        writeBin(1, outfile)
        stop("Interrupted normalization")
    }, .package = "moleculaR")
    input$normalization <- "TIC"
    expect_error(.shinyInitialize(state, input, "upload", NULL), "Interrupted normalization")
    partial <- state$normalizedPath
    expect_true(file.exists(partial))
    input$normalization <- "none"
    .shinyInitialize(state, input, "upload", NULL)
    expect_false(file.exists(partial))
})

test_that("Shiny rejects unsupported MSI and validates run, ROI and profile model separately", {
    state <- shiny_state()
    on.exit(.shinyRemove(state$root, state$root, TRUE))
    for (status in c(FALSE, NA)) {
        x <- make_msi_arrays(centroided = status)
        path <- write_msi_pair(x, tempfile(fileext = ".imzML"))
        expect_error(.shinyInitialize(state, shiny_input(shiny_upload(path)), "upload", NULL), "tagged centroided")
        expect_false(state$ready)
    }
    x <- make_msi_arrays(runs = c("a", "b", "a"))
    # imzML encodes one acquisition; use the imported native multi-run object to
    # verify the selection without pretending its writer retains run labels.
    path <- write_msi_pair(x, tempfile(fileext = ".imzML"))
    testthat::local_mocked_bindings(readMsi = function(...) x, .package = "moleculaR")
    input <- shiny_input(shiny_upload(path))
    expect_error(.shinyInitialize(state, input, "upload", NULL), "Select one experimental run")
    input$run <- "a"
    .shinyInitialize(state, input, "upload", NULL)
    expect_equal(.msiNumberOfPixels(state$msData), 2)
    expect_equal(unique(.msiRun(state$msData)), "a")
    expect_equal(spatstat.geom::area.owin(state$spwin), 2)
    expect_length(Cardinal::mz(state$msData)[[2]], 0)
    input$useRoi <- TRUE
    input$xMin <- input$xMax <- 3
    input$yMin <- input$yMax <- 1
    .shinyInitialize(state, input, "upload", NULL)
    expect_equal(.msiNumberOfPixels(state$msData), 1)
    expect_equal(spatstat.geom::area.owin(state$spwin), 1)
    input$xMin <- input$xMax <- 999
    expect_error(.shinyInitialize(state, input, "upload", NULL), "no acquired pixels")
    input$xMin <- 3
    input$xMax <- 2
    expect_error(.shinyInitialize(state, input, "upload", NULL), "valid x/y ROI")
    settings <- .shinyDataSettings(shiny_input(fwhmSource = "spectrum"))
    expect_error(.shinyModel(settings, NULL), "separate profile")
    profile <- system.file("extdata", "pos-XIII-82492-singleSpectrum.tsv", package = "moleculaR")
    measured <- .shinyModel(settings, profile)
    expect_s3_class(measured$model, "fwhm")
    expect_true(all(is.finite(getFwhm(measured$model, c(400, 800)))))
    expect_match(measured$description, "separate uploaded profile")
    expect_identical(measured$model$estimation$peakSampling, .Machine$integer.max)
    for (modality in c("FTICR", "TIMS-TOF", "TOF-reflector")) {
        model <- .shinyModel(.shinyDataSettings(shiny_input(fwhmModality = modality, customMassRes = FALSE)), NULL)
        expect_equal(getFwhm(model$model, 249), 249 / switch(modality, FTICR = 123000, `TIMS-TOF` = 40000, `TOF-reflector` = 10000))
    }
})

test_that("Shiny maps match direct native APIs and downloads remain portable", {
    state <- shiny_state()
    on.exit(.shinyRemove(state$root, state$root, TRUE))
    path <- write_msi_pair(fixture_cardinal(), tempfile(fileext = ".imzML"))
    input <- shiny_input(shiny_upload(path))
    .shinyInitialize(state, input, "upload", NULL)
    result <- suppressWarnings(.shinyAnalyze(state, input, "single"))
    expected <- batchAnalyteSearch(state$msData, state$fwhmObj, state$spwin, input$mz, verbose = FALSE)
    expect_equal(shiny_fields(result$pointPattern), shiny_fields(expected), ignore_attr = TRUE)
    expect_equal(result$pointPattern$metaData$mzVals, input$mz)
    reference <- suppressWarnings(probMap(expected, bw = 1, seed = 1, verbose = FALSE))
    for (field in c("rhoMoi", "rhoCsr", "hotspotIm", "coldspotIm")) expect_equal(result$map[[field]]$v, reference[[field]]$v)
    legacyInput <- input
    legacyInput$mz <- 100
    historical <- suppressWarnings(.shinyAnalyze(state, legacyInput, "single"))
    for (field in c("rhoMoi", "rhoCsr")) expect_equal(historical$map[[field]]$v, legacy_fixture()$mpm[[field]], tolerance = 1e-10)
    expect_equal(historical$map$hotspotIm$v, legacy_fixture()$mpm$hotspot, tolerance = 1e-10)
    expect_equal(historical$map$coldspotIm$v, legacy_fixture()$mpm$coldspot, tolerance = 1e-10)
    batch <- suppressWarnings(.shinyAnalyze(state, input, "batch"))
    expected <- batchAnalyteSearch(state$msData, state$fwhmObj, state$spwin, c(100, 99.9), verbose = FALSE)
    expect_equal(shiny_fields(batch$pointPattern), shiny_fields(expected), ignore_attr = TRUE)
    expect_equal(nrow(batch$pointPattern$metaData), 2)
    input$standardize <- TRUE
    standardized <- suppressWarnings(.shinyAnalyze(state, input, "batch"))
    expect_equal(shiny_fields(standardized$pointPattern), shiny_fields(transformIntensity(expected, "z-score")), ignore_attr = TRUE)
    for (format in c("rds", "csv", "png")) {
        file <- tempfile(fileext = paste0(".", format))
        .shinySaveAnalysis(result, file, format)
        expect_gt(file.size(file), 100)
        if (format == "rds") {
            restored <- readRDS(file)
            expect_equal(restored$map$rhoMoi$v, result$map$rhoMoi$v)
            expect_false(any(grepl(state$root, capture.output(str(restored)), fixed = TRUE)))
            expect_null(restored$msData)
        } else if (format == "csv") {
            table <- utils::read.csv(file)
            expect_equal(table$intensity, result$pointPattern$marks$intensity)
            expect_equal(table$mzVals, rep(input$mz, nrow(table)))
        } else expect_identical(readBin(file, "raw", 8), as.raw(c(137, 80, 78, 71, 13, 10, 26, 10)))
    }
    plotFile <- tempfile(fileext = ".pdf")
    grDevices::pdf(plotFile, width = 12, height = 12)
    for (view in c("MPM", "image", "detailed")) expect_no_error(.shinyPlotAnalysis(result, view))
    grDevices::dev.off()
    expect_gt(file.size(plotFile), 100)
    expect_error(.shinySaveAnalysis(NULL, tempfile()), "no current analysis")
    input$mz <- 10000
    expect_error(.shinyAnalyze(state, input, "single"), "mass range")
    input$mz <- c(100, 101)
    expect_error(.shinyAnalyze(state, input, "single"), "one theoretical")
    input <- shiny_input(confirmedOnly = TRUE)
    expect_error(.shinyAnalyze(state, input, "single"), "Upload verified masses")
    state$verifiedMasses <- numeric()
    expect_error(.shinyAnalyze(state, input, "single"), "No analyte signal")
    state$ready <- FALSE
    expect_error(.shinyAnalyze(state, input, "single"), "initialize")
})

test_that("Shiny lipid groups preserve annotations and explicit standardization errors", {
    database <- data.frame("Lipid ID" = c("a", "b", "c", "d"), "Formula (pH7.3)" = "C10H20",
        "Abbreviation*" = c("PC(30:0)", "PC(32:1)", "PE(34:2)", "PE(36:3)"), numDoubleBond = 0:3,
        lipidGroup = c("PC(x:x)", "PC(x:x)", "PE(x:x)", "PE(x:x)"), chainLength = 30:33,
        "Exact m/z of [M-H]-" = 100, "Exact m/z of [M+H]+" = 100.1,
        "Exact m/z of [M+Na]+" = 99.9, "Exact m/z of [M+K]+" = 100, check.names = FALSE)
    state <- shiny_state()
    on.exit(.shinyRemove(state$root, state$root, TRUE))
    path <- write_msi_pair(fixture_cardinal(), tempfile(fileext = ".imzML"))
    input <- shiny_input(shiny_upload(path))
    .shinyInitialize(state, input, "upload", NULL)
    state$lipids <- database
    for (mode in c("class", "ion", "saturation")) {
        input$lipidView <- mode
        result <- suppressWarnings(.shinyAnalyze(state, input, "lipid"))
        expect_equal(nrow(result$pointPattern$metaData), switch(mode, class = 2L, ion = 4L, saturation = 1L))
        expect_true(all(result$pointPattern$metaData$adduct == "M-H"))
        expect_equal(result$pointPattern$metaData$lipidID, switch(mode, class = c("a", "b"), ion = c("a", "b", "c", "d"), saturation = "b"))
    }
    for (group in c("sat", "mono", "di", "poly")) expect_equal(nrow(.shinyLipidSelection(database, "saturation", NULL, NULL, group, "M-H")$database), 1)
    input$lipidView <- "class"
    input$lipidClass <- "absent"
    expect_error(.shinyAnalyze(state, input, "lipid"), "Select a lipid group")
    constant <- searchAnalyte(100, 0.1, state$msData, spwin = state$spwin)
    constant$marks$intensity[] <- 1
    expect_error(.shinyStandardize(constant), "constant or insufficient")
})

test_that("Shiny session actions invalidate correctly and close only owned files", {
    path <- write_msi_pair(fixture_cardinal(), tempfile(fileext = ".imzML"))
    upload <- shiny_upload(path)
    ownRoot <- NULL
    shiny::testServer(.cardinalAnalysisServer, args = list(source = "upload"), {
        session$flushReact()
        do.call(session$setInputs, shiny_input(upload))
        session$setInputs(loadData = 1)
        expect_true(state$ready)
        ownRoot <<- state$root
        ownPath <- state$sourcePath
        session$setInputs(renderMPM = 1)
        expect_false(is.null(state$result))
        expect_match(output$message, "complete")
        map <- state$result$map$rhoMoi$v
        cache <- state$cache
        session$setInputs(plotView = "image")
        expect_equal(state$result$map$rhoMoi$v, map)
        expect_identical(state$cache, cache)
        session$setInputs(mz = 100.001)
        expect_null(state$result)
        expect_true(state$ready)
        session$setInputs(renderMPM = 2)
        expect_false(is.null(state$result))
        session$setInputs(fwhmMassRes = 2000)
        expect_false(state$ready)
        expect_null(state$result)
        expect_length(ls(state$cache), 0)
        session$setInputs(renderMPM = 3)
        expect_match(output$message, "initialize")
        session$setInputs(loadData = 2)
        expect_true(state$ready)
        expect_equal(getFwhm(state$fwhmObj, 100), 0.05)
        session$setInputs(loadLipids = 1)
        expect_equal(nrow(state$lipids), 10319)
        database <- state$lipids[1:2, setdiff(names(state$lipids), c("numDoubleBond", "lipidGroup", "chainLength")), drop = FALSE]
        database$Level <- "Species"
        lipidPath <- tempfile(fileext = ".tsv")
        utils::write.table(database, lipidPath, sep = "\t", row.names = FALSE, quote = FALSE)
        session$setInputs(lipidFile = data.frame(name = "species.tsv", datapath = lipidPath))
        expect_null(state$lipids)
        session$setInputs(loadLipids = 2)
        expect_equal(nrow(state$lipids), 2)
        expect_true(all(.shinyLipidColumns() %in% names(state$lipids)))
        expect_true(file.exists(lipidPath))
        # A concurrent session cannot use this model, results, or owned files.
        other <- shiny::MockShinySession$new()
        shiny::withReactiveDomain(other, {
            otherModule <- .cardinalAnalysisServer("other", "upload")
            other$flushReact()
            expect_false(otherModule$state$ready)
            expect_false(identical(otherModule$state$root, ownRoot))
            inputs <- shiny_input(upload)
            names(inputs) <- paste0("other-", names(inputs))
            do.call(other$setInputs, inputs)
            other$setInputs(`other-renderMPM` = 1)
            expect_match(otherModule$state$message, "initialize")
            other$setInputs(`other-loadData` = 1)
            expect_true(otherModule$state$ready)
            expect_true(file.exists(ownPath))
            other$setInputs(`other-resetData` = 1)
            expect_false(otherModule$state$ready)
            expect_null(otherModule$state$msData)
            expect_true(file.exists(ownPath))
            other$close()
        })
        expect_true(file.exists(ownPath))
        expect_true(all(file.exists(upload$datapath)))
    })
    expect_false(dir.exists(ownRoot))
    expect_true(all(file.exists(upload$datapath)))
})

test_that("both Shiny interfaces use namespaced native workflows and validate launch choices", {
    for (only in c(FALSE, TRUE)) {
        html <- htmltools::renderTags(.cardinalAppUi(only))$html
        expect_match(html, "Separate profile spectrum")
        expect_match(html, "Aggregated intensity image")
        expect_match(html, "example-loadData")
        expect_identical(grepl("main-loadData", html, fixed = TRUE), !only)
        ids <- regmatches(html, gregexpr(' id="[^"]+"', html))[[1]]
        expect_false(anyDuplicated(ids) > 0)
    }
    for (bad in list(NULL, NA_character_, "www", "missing", c("package-app", "example-app"), 1)) expect_error(runShiny(bad), "valid app")
    expect_error(runShiny(), "valid app")
})
