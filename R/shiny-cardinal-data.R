# Native Shiny data operations. Session paths and numeric metadata stay separate
# from the Cardinal experiment; no complete spectral matrix or legacy list is built.

.shinyValue <- function(value, default) if (is.null(value) || !length(value)) default else value
.shinyCache <- function() new.env(parent = emptyenv())

.shinySessionDirectory <- function() {
    root <- tempfile("moleculaR-session-", tmpdir = tempdir())
    dir.create(root)
    file.create(file.path(root, ".moleculaR-session"))
    normalizePath(root, winslash = "/", mustWork = TRUE)
}

.shinyOwnedPath <- function(path, root, allowRoot = FALSE) {
    root <- normalizePath(root, winslash = "/", mustWork = TRUE)
    parent <- normalizePath(tempdir(), winslash = "/", mustWork = TRUE)
    path <- normalizePath(path, winslash = "/", mustWork = FALSE)
    if (!identical(dirname(root), parent) || !grepl("^moleculaR-session-", basename(root)) ||
        !file.exists(file.path(root, ".moleculaR-session")) ||
        !(startsWith(path, paste0(root, "/")) || (allowRoot && identical(path, root)))) {
        stop("Refusing cleanup outside the owned Shiny session directory.")
    }
    path
}

.shinyRemove <- function(path, root, allowRoot = FALSE) {
    if (is.null(path) || !length(path) || !file.exists(path)) return(invisible(NULL))
    target <- .shinyOwnedPath(path, root, allowRoot)
    unlink(target, recursive = isTRUE(file.info(target)$isdir))
    invisible(NULL)
}

# Stamps identify Shiny uploads, not the continued existence of their temporary
# paths. Once copied, a session can keep using a file after upload paths expire.
.shinyUploadStamp <- function(file) {
    if (is.null(file)) return(NULL)
    as.list(file[, intersect(c("name", "datapath", "size"), names(file)), drop = FALSE])
}

.shinyCopyPair <- function(upload, root) {
    if (!is.data.frame(upload) || nrow(upload) != 2L ||
        !all(c("name", "datapath") %in% names(upload))) {
        stop("Upload one matching imzML and ibd pair.")
    }
    extensions <- tolower(tools::file_ext(upload$name))
    if (!setequal(extensions, c("imzml", "ibd")) || anyDuplicated(extensions)) {
        stop("Upload exactly one imzML file and one ibd file.")
    }
    stems <- tolower(tools::file_path_sans_ext(basename(upload$name)))
    if (length(unique(stems)) != 1L) stop("The imzML and ibd filenames must have the same stem.")
    if (any(!file.exists(upload$datapath)) || any(file.info(upload$datapath)$isdir)) {
        stop("An uploaded file is unavailable; upload the pair again.")
    }
    directory <- tempfile("dataset-", tmpdir = root)
    dir.create(directory)
    paths <- file.path(directory, c("dataset.imzML", "dataset.ibd"))
    copied <- file.copy(upload$datapath[match(c("imzml", "ibd"), extensions)], paths)
    if (!all(copied)) {
        .shinyRemove(directory, root)
        stop("Could not prepare the uploaded imzML/ibd pair.")
    }
    list(path = paths[1], directory = directory,
         filename = basename(upload$name[extensions == "imzml"]))
}

.shinyAuxiliary <- function(state, upload, role, extension) {
    if (is.null(upload)) return(NULL)
    if (!is.data.frame(upload) || nrow(upload) != 1L ||
        !all(c("name", "datapath") %in% names(upload)) ||
        tolower(tools::file_ext(upload$name)) != extension) {
        stop("Upload one ", extension, " file for ", role, ".")
    }
    stamp <- .shinyUploadStamp(upload)
    previous <- state$auxiliary[[role]]
    if (!is.null(previous) && identical(previous$stamp, stamp) && file.exists(previous$path)) return(previous$path)
    if (!file.exists(upload$datapath)) stop("The uploaded ", role, " is unavailable; upload it again.")
    path <- tempfile(paste0(role, "-"), tmpdir = state$root, fileext = paste0(".", extension))
    if (!file.copy(upload$datapath, path)) stop("Could not copy the uploaded ", role, ".")
    if (!is.null(previous)) .shinyRemove(previous$path, state$root)
    auxiliary <- state$auxiliary
    auxiliary[[role]] <- list(path = path, stamp = stamp)
    state$auxiliary <- auxiliary
    path
}

.shinyNumbers <- function(text, label = "m/z values") {
    if (!is.character(text) || length(text) != 1L || is.na(text) || !nzchar(trimws(text))) stop("Enter ", label, ".")
    if (grepl("(^[,;]|[,;][[:space:]]*[,;]|[,;][[:space:]]*$)", trimws(text))) stop("Missing entry in ", label, ".")
    tokens <- strsplit(trimws(text), "[,;[:space:]]+")[[1]]
    numbers <- suppressWarnings(as.numeric(tokens))
    if (!length(numbers) || any(!is.finite(numbers)) || any(numbers <= 0)) stop(label, " must be positive finite numbers.")
    numbers
}

.shinyAnnotations <- function(path) {
    # Ordinary mz CSVs and METASPACE exports with two metadata rows are accepted.
    for (skip in c(0L, 2L)) {
        table <- tryCatch(utils::read.csv(path, skip = skip, check.names = FALSE), error = function(e) NULL)
        if (is.null(table)) next
        column <- which(tolower(trimws(names(table))) == "mz")
        if (length(column) == 1L) {
            masses <- suppressWarnings(as.numeric(as.character(table[[column]])))
            .msiValidateMasses(masses)
            return(masses)
        }
    }
    stop("Annotations need a CSV column named mz, pre-filtered to your chosen FDR.")
}

.shinyMzRange <- function(mz) {
    .msiValidateMasses(mz)
    if (length(mz)) range(mz) else c(Inf, -Inf)
}

.shinyMassLimits <- function(x) {
    if (methods::is(x, "MSImagingExperiment")) {
        limits <- .shinyMzRange(Cardinal::mz(x))
    } else {
        # Read only the m/z arrays; never unlist all masses or read intensities.
        bounds <- Cardinal::spectrapply(x, .shinyMzRange, spectra = "mz", simplify = TRUE,
            verbose = FALSE, BPPARAM = BiocParallel::SerialParam())
        bounds <- matrix(as.numeric(bounds), nrow = 2)
        limits <- c(min(bounds[1, ]), max(bounds[2, ]))
    }
    if (any(!is.finite(limits))) stop("The selected data contain no valid measured masses.")
    limits
}

.shinyDataSettings <- function(input, runs = character()) {
    source <- .shinyValue(input$fwhmSource, "default")
    run <- .shinyValue(input$run, "")
    if (length(runs) == 1L) run <- runs[1]
    roi <- if (isTRUE(input$useRoi)) c(input$xMin, input$xMax, input$yMin, input$yMax) else NULL
    spacing <- trimws(.shinyValue(input$pixelSpacing, ""))
    list(run = run, roi = roi, useRoi = isTRUE(input$useRoi), pixelSpacing = spacing,
        checkBinary = isTRUE(input$checkBinary), normalization = .shinyValue(input$normalization, "none"),
        scale = if (identical(input$normalizationScale, "unit")) 1 else NA_real_,
        standardMz = if (identical(input$normalization, "IS")) input$standardMz else NULL,
        standardTol = if (identical(input$normalization, "IS")) input$standardTol else NULL,
        fwhmSource = source,
        modality = if (source == "default") .shinyValue(input$fwhmModality, "FTICR") else NULL,
        resolution = if (source == "default" && isTRUE(input$customMassRes)) c(massRes = input$fwhmMassRes, mass = input$fwhmMass) else NULL,
        header = if (source == "spectrum") isTRUE(input$profileHeader) else NULL,
        peakMethod = if (source == "spectrum") .shinyValue(input$peakMethod, "diff") else NULL,
        SNR = if (source == "spectrum") .shinyValue(input$profileSnr, 5) else NULL,
        profile = if (source == "spectrum") .shinyUploadStamp(input$spectrFile) else NULL,
        annotations = .shinyUploadStamp(input$annotationFile))
}

.shinyModel <- function(settings, profilePath) {
    if (identical(settings$fwhmSource, "default")) {
        if (!is.null(settings$resolution) && length(settings$resolution) != 2L) stop("Supply resolving power and reference m/z.")
        model <- suppressMessages(estimateFwhm(modality = settings$modality, userMassRes = settings$resolution))
        reference <- if (is.null(settings$resolution)) 249 else unname(settings$resolution["mass"])
        description <- paste0("Approximate ", settings$modality, "; resolving power ",
            signif(reference / getFwhm(model, reference), 7), " at m/z ", reference, ".")
    } else if (identical(settings$fwhmSource, "spectrum")) {
        if (is.null(profilePath)) stop("Upload a separate profile spectrum TSV to estimate measured widths.")
        model <- suppressMessages(estimateFwhm(readSingleSpect(profilePath, header = settings$header),
            peakMethod = settings$peakMethod, SNR = settings$SNR, peakSampling = .Machine$integer.max))
        description <- paste0("FWHM estimated from the separate uploaded profile spectrum; Cardinal ", settings$peakMethod, ".")
    } else stop("Choose an instrument model or a separate profile spectrum.")
    list(model = model, description = description)
}

.shinyInitialize <- function(state, input, source, examplePath, progress = function(...) NULL) {
    state$ready <- FALSE
    state$result <- NULL
    state$cache <- .shinyCache()
    state$msData <- NULL
    .shinyRemove(state$normalizedPath, state$root)
    state$normalizedPath <- NULL
    progress(0.1, "Preparing files")
    if (source == "upload") {
        stamp <- .shinyUploadStamp(input$imzmlFile)
        if (is.null(state$sourcePath) || !identical(stamp, state$uploadStamp)) {
            pair <- .shinyCopyPair(input$imzmlFile, state$root)
            .shinyRemove(state$datasetDirectory, state$root)
            state$sourcePath <- pair$path
            state$datasetDirectory <- pair$directory
            state$sourceName <- pair$filename
            state$uploadStamp <- stamp
        }
    } else {
        state$sourcePath <- file.path(examplePath, "example.imzML")
        state$sourceName <- "Reduced glioblastoma example (32 x 32 coordinate crop)"
        if (is.null(input$lipidFile)) state$lipids <- readRDS(file.path(examplePath, "lipids.rds"))
    }
    settings <- .shinyDataSettings(input)
    progress(0.2, "Reading centroided MSI")
    candidate <- readMsi(state$sourcePath, memory = FALSE, as = "MSImagingArrays",
                         check = settings$checkBinary, verbose = FALSE)
    state$runs <- unique(.msiRun(candidate))
    if (length(state$runs) == 1L) settings$run <- state$runs[1]
    if (!nzchar(settings$run) || !settings$run %in% state$runs) stop("Select one experimental run, then initialize again.")
    candidate <- .msiSelectPixels(candidate, which(.msiRun(candidate) == settings$run))
    if (!isTRUE(ProtGenerics::centroided(candidate))) {
        stop("MSI must be tagged centroided. Peak-pick profile MSI explicitly in Cardinal before upload; a separate profile TSV only estimates widths.")
    }
    coordinates <- .msiCoordinates(candidate)
    state$coordinateBounds <- lapply(coordinates, range)
    if (settings$useRoi) {
        bounds <- settings$roi
        if (!is.numeric(bounds) || length(bounds) != 4L || any(!is.finite(bounds)) || bounds[1] > bounds[2] || bounds[3] > bounds[4]) stop("Enter valid x/y ROI bounds.")
        pixels <- which(coordinates$x >= bounds[1] & coordinates$x <= bounds[2] & coordinates$y >= bounds[3] & coordinates$y <= bounds[4])
        if (!length(pixels)) stop("The ROI contains no acquired pixels.")
        candidate <- .msiSelectPixels(candidate, pixels)
    }
    .msiValidate(candidate, extraction = TRUE)
    spacing <- if (nzchar(settings$pixelSpacing)) .shinyNumbers(settings$pixelSpacing, "pixel spacing") else NULL
    if (!is.null(spacing) && !length(spacing) %in% c(1L, 2L)) stop("Pixel spacing needs one or two values.")
    progress(0.35, "Preparing tissue window and FWHM model")
    window <- tryCatch(createSpatialWindow(.msiCoordinates(candidate), pixelSize = spacing, plot = FALSE),
        error = function(e) stop("Cannot infer the tissue grid. Supply pixel spacing or choose a larger ROI. ", conditionMessage(e)))
    profile <- if (settings$fwhmSource == "spectrum") .shinyAuxiliary(state, input$spectrFile, "profile", "tsv") else NULL
    model <- .shinyModel(settings, profile)
    annotationPath <- .shinyAuxiliary(state, input$annotationFile, "annotations", "csv")
    annotations <- if (!is.null(annotationPath)) .shinyAnnotations(annotationPath) else if (source == "example") .shinyAnnotations(file.path(examplePath, "annotations.csv")) else NA_real_
    progress(0.5, "Inspecting measured mass range")
    limits <- .shinyMassLimits(candidate)
    if (settings$normalization != "none") {
        progress(0.65, "Normalizing intensities")
        path <- tempfile("normalized-", tmpdir = state$root, fileext = ".bin")
        # Track the owned target before writing so a failed/partial write is
        # retired on the next initialization as well as on session closure.
        state$normalizedPath <- path
        candidate <- normalizeIntensity(candidate, settings$normalization, scale = settings$scale,
            mz = settings$standardMz, tol = settings$standardTol, outfile = path)
    }
    state$msData <- candidate
    state$fwhmObj <- model$model
    state$modelDescription <- model$description
    state$spwin <- window
    state$massLimits <- limits
    state$verifiedMasses <- annotations
    state$appliedSettings <- settings
    state$ready <- TRUE
    state$message <- paste0("Ready: ", .msiNumberOfPixels(candidate), " acquired pixels; m/z ",
        paste(signif(limits, 7), collapse = " to "), ". ", model$description)
    progress(1, "Ready")
    invisible(TRUE)
}

.shinyLipidColumns <- function() c("Lipid ID", "Formula (pH7.3)", "Abbreviation*", "numDoubleBond", "lipidGroup", "chainLength",
    "Exact m/z of [M-H]-", "Exact m/z of [M+H]+", "Exact m/z of [M+Na]+", "Exact m/z of [M+K]+")

.shinyLipidSelection <- function(database, mode, lipidClass, ion, saturation, adducts) {
    if (!all(.shinyLipidColumns() %in% names(database))) stop("The lipid database lacks required species metadata.")
    if (mode == "class") {
        selected <- database$lipidGroup == lipidClass
    } else {
        groups <- c("LPA(x:x)", "LPC(x:x)", "LPE(x:x)", "LPG(x:x)", "LPI(x:x)", "LPS(x:x)",
                    "PA(x:x)", "PC(x:x)", "PE(x:x)", "PG(x:x)", "PI(x:x)", "PS(x:x)")
        selected <- database$lipidGroup %in% groups
        if (mode == "ion") adducts <- ion else if (mode == "saturation") {
            count <- database$numDoubleBond
            selected <- selected & switch(saturation, sat = count == 0, mono = count == 1,
                di = count == 2, poly = count > 2, stop("Choose a saturation group."))
        } else stop("Choose a lipid grouping.")
    }
    labels <- c("[M-H]-" = "M-H", "[M+H]+" = "M+H", "[M+Na]+" = "M+Na", "[M+K]+" = "M+K")
    adducts <- unname(ifelse(adducts %in% names(labels), labels[adducts], adducts))
    list(database = database[which(selected), , drop = FALSE], adducts = adducts)
}

.shinyStandardize <- function(pattern) {
    if (!pattern$n) return(pattern)
    for (id in pattern$metaData$idx) {
        values <- pattern$marks$intensity[pattern$marks$idx == id]
        if (length(values) < 2L || !is.finite(stats::sd(values)) || stats::sd(values) <= 0) {
            stop("Some analytes have constant or insufficient signal. Disable standardization or remove those analytes.")
        }
    }
    transformIntensity(pattern, method = "z-score")
}

.shinyAnalyze <- function(state, input, kind, progress = function(...) NULL) {
    if (!isTRUE(state$ready) || is.null(state$msData)) stop("Load and initialize the data and FWHM model first.")
    method <- .shinyValue(input$wMethod, "Gaussian")
    confirmed <- isTRUE(input$confirmedOnly)
    if (confirmed && length(state$verifiedMasses) == 1L && is.na(state$verifiedMasses)) stop("Upload verified masses before selecting confirmed-only results.")
    progress(0.1, "Extracting analyte intensities")
    if (kind == "single") {
        masses <- .shinyValue(input$mz, 788.5447)
        if (length(masses) != 1L) stop("Enter one theoretical m/z value.")
    } else if (kind == "batch") masses <- unique(.shinyNumbers(input$mzStr)) else masses <- NULL
    if (!is.null(masses)) {
        if (!is.numeric(masses) || any(!is.finite(masses)) || any(masses <= 0) || any(masses < state$massLimits[1] | masses > state$massLimits[2])) stop("Theoretical m/z values must be within the selected data's mass range.")
        pattern <- batchAnalyteSearch(state$msData, state$fwhmObj, spwin = state$spwin, m = masses,
            wMethod = method, verifiedMasses = state$verifiedMasses, confirmedOnly = confirmed,
            cache = state$cache, verbose = FALSE)
    } else if (kind == "lipid") {
        if (is.null(state$lipids)) stop("Load a lipid database first.")
        selection <- .shinyLipidSelection(state$lipids, .shinyValue(input$lipidView, "class"),
            input$lipidClass, input$lipidIon, .shinyValue(input$lipidSat, "sat"), input$adducts)
        if (!nrow(selection$database) || !length(selection$adducts)) stop("Select a lipid group and at least one adduct.")
        pattern <- batchLipidSearch(state$msData, state$fwhmObj, spwin = state$spwin,
            sldb = selection$database, adduct = selection$adducts, wMethod = method,
            verifiedMasses = state$verifiedMasses, confirmedOnly = confirmed, cache = state$cache, verbose = FALSE)
    } else stop("Unknown analysis type.")
    if (!pattern$n) stop("No analyte signal was detected for this query and confirmation setting.")
    if (kind != "single" && isTRUE(input$standardize)) pattern <- .shinyStandardize(pattern)
    seed <- .shinyValue(input$seed, 1)
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) || seed < 0 || seed != floor(seed) || seed > .Machine$integer.max) stop("Seed must be a non-negative integer.")
    bandwidth <- if (identical(input$bandwidthMode, "fixed")) input$bandwidth else seq(0.5, 10)
    if (!is.numeric(bandwidth) || !length(bandwidth) || any(!is.finite(bandwidth)) || any(bandwidth <= 0)) stop("Supply a positive bandwidth.")
    progress(0.55, "Calculating the probability map")
    map <- probMap(pattern, bw = bandwidth, seed = seed, verbose = FALSE)
    progress(0.9, "Preparing image and downloads")
    settings <- state$appliedSettings
    # Export only portable scientific settings, not session/upload file paths.
    settings$profile <- settings$annotations <- NULL
    list(map = map, pointPattern = pattern, ionImage = spp2im(pattern, rescale = FALSE),
        fwhm = state$fwhmObj,
        provenance = list(source = state$sourceName, kind = kind, settings = settings,
            requestedMasses = masses, aggregation = method, confirmedOnly = confirmed,
            lipidSelection = if (kind == "lipid") list(mode = input$lipidView, class = input$lipidClass,
                ion = input$lipidIon, saturation = input$lipidSat, adducts = selection$adducts) else NULL,
            verifiedMasses = state$verifiedMasses, modelDescription = state$modelDescription,
            versions = c(moleculaR = as.character(utils::packageVersion("moleculaR")),
                Cardinal = as.character(utils::packageVersion("Cardinal"))),
            standardized = kind != "single" && isTRUE(input$standardize), seed = seed,
            bandwidth = unname(map$bw), acquiredPixels = .msiNumberOfPixels(state$msData)))
}

.shinyPointTable <- function(result) {
    if (is.null(result)) stop("There is no current analysis to download.")
    pattern <- result$pointPattern
    metadata <- pattern$metaData[match(pattern$marks$idx, pattern$metaData$idx), , drop = FALSE]
    rownames(metadata) <- NULL
    data.frame(x = pattern$x, y = pattern$y, intensity = pattern$marks$intensity,
               metadata, check.names = FALSE)
}

.shinySaveAnalysis <- function(result, file, format = c("rds", "csv", "png")) {
    if (is.null(result)) stop("There is no current analysis to download.")
    format <- match.arg(format)
    if (format == "rds") saveRDS(result, file) else if (format == "csv") {
        utils::write.csv(.shinyPointTable(result), file, row.names = FALSE)
    } else {
        grDevices::png(file, width = 1200, height = 1000, res = 130)
        on.exit(grDevices::dev.off(), add = TRUE)
        graphics::plot(result$map, what = "MPM")
    }
    invisible(file)
}
