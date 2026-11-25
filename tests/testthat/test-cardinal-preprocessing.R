phase4_fixture <- function() readRDS(test_path("..", "fixtures", "legacy", "phase4", "preprocessing.rds"))

phase4_input <- function(arrays = TRUE, dropZeros = TRUE) {
    f <- phase4_fixture()
    if (!arrays) return(make_msi_experiment(f$axis, f$signal, f$coordinates))
    keep <- lapply(seq_len(ncol(f$signal)), function(i) if (dropZeros) f$signal[, i] > 0 else rep(TRUE, length(f$axis)))
    make_msi_arrays(lapply(keep, function(k) f$axis[k]),
        lapply(seq_along(keep), function(i) f$signal[keep[[i]], i]), f$coordinates)
}

phase4_values <- function(x) {
    if (methods::is(x, "MSImagingArrays")) {
        lapply(seq_len(.msiNumberOfPixels(x)), function(i) as.numeric(Cardinal::intensity(x)[[i]]))
    } else as.matrix(Cardinal::intensity(x))
}

test_that("FTICR calibration is numeric and preserves knots, defaults and extrapolation", {
    f <- phase4_fixture()
    reference <- readRDS(system.file("extdata", "fticr-fwhm-reference.rds", package = "moleculaR"))
    expect_true(all(vapply(reference$interpolation, is.atomic, logical(1))))
    expect_equal(getFwhm(useDefaultFwhm(), f$fticrMasses), f$fticrPredictions, tolerance = 1e-12)
    legacy <- new.env(parent = asNamespace("moleculaR"))
    sys.source(test_path("..", "fixtures", "legacy", "phase4", "useDefaultFwhm.R"), legacy)
    grid <- sort(unique(c(seq(50, 2000, length.out = 2001), reference$interpolation$x)))
    for (anchor in list(NULL, c(massRes = 150000, mass = 400))) {
        expect_equal(getFwhm(useDefaultFwhm(anchor), grid),
                     getFwhm(legacy$useDefaultFwhm(anchor), grid), tolerance = 1e-12)
    }
    # A default model must not load any part of the legacy dataset.
    local_mocked_bindings(data = function(...) stop("legacy data load"), .package = "utils")
    expect_s3_class(useDefaultFwhm(), "fwhm")
})

test_that("single-spectrum TSV input is numeric, validated and independent of MALDIquant", {
    path <- tempfile(fileext = ".tsv")
    writeLines(c("mz\tintensity", "100\t0", "101\t2", "102\t1"), path)
    x <- readSingleSpect(path)
    expect_s3_class(x, "moleculaRSpectrum")
    expect_equal(x$mz, c(100, 101, 102))
    expect_equal(x$intensity, c(0, 2, 1))
    writeLines(c("100\t0", "101\t2"), path)
    expect_equal(readSingleSpect(path, header = FALSE)$mz, c(100, 101))
    for (lines in list(c("mz\ti", "100\tNA"), c("mz\ti", "100\tInf"),
        c("mz\ti", "101\t1", "100\t2"), c("mz\ti", "100\t1", "100\t2"),
        c("mz\ti", "0\t1"), c("mz\ti", "100\t-1"), c("mz\ti", "bad\t2"),
        c("mz\ti\tthird", "100\t2\t0"), "mz\ti")) {
        writeLines(lines, path)
        expect_error(readSingleSpect(path))
    }
    expect_error(readSingleSpect("absent-spectrum.tsv"), "path")
    expect_error(readSingleSpect(path, header = NA), "header")
})

test_that("numeric half-height widths reproduce frozen legacy calculations", {
    f <- phase4_fixture()
    widths <- .estimatePeakWidths(f$profile$mz, f$profile$intensity, f$peakIndices)
    expect_equal(widths, f$widths$fwhmValues, tolerance = 1e-12)
    for (input in list(f$profile, as.data.frame(f$profile),
                       .numericSpectrum(f$profile$mz, f$profile$intensity))) {
        model <- estimateFwhm(input, peakIndices = f$peakIndices)
        expect_equal(getFwhm(model, f$predictionMasses), f$predictions, tolerance = 1e-10)
        expect_identical(model$estimation$peakSelection, "supplied peakIndices")
    }
    model <- estimateFwhm(mz = f$profile$mz, intensity = f$profile$intensity, peakIndices = f$peakIndices)
    expect_equal(getFwhm(model, f$predictionMasses), f$predictions, tolerance = 1e-10)
    expect_true(all(is.na(.estimatePeakWidths(1:3, c(2, 3, 2), c(1, 2, 3)))))
    expect_error(.estimatePeakWidths(1:3, c(0, 1, 0), c(2, 2)), "peakIndices")
    expect_error(.estimatePeakWidths(1:3, c(0, 1, 0), 4), "peakIndices")
    expect_error(estimateFwhm(f$profile, peakIndices = f$peakIndices[1:4]), "five distinct")
})

test_that("Cardinal profile picking and native precomputed widths support both representations", {
    f <- phase4_fixture()
    expect_message(model <- estimateFwhm(f$profile), "Cardinal diff")
    expect_equal(length(model$peaks), 21L)
    expect_equal(getFwhm(model, f$predictionMasses), f$predictions, tolerance = 1e-8)
    for (arrays in c(TRUE, FALSE)) {
        profile <- if (arrays) make_msi_arrays(rep(list(f$profile$mz), 2), rep(list(f$profile$intensity), 2),
            data.frame(x = 1:2, y = 1), centroided = FALSE) else
            Cardinal::MSImagingExperiment(spectraData = list(intensity = cbind(f$profile$intensity, f$profile$intensity)),
                featureData = Cardinal::MassDataFrame(f$profile$mz),
                pixelData = Cardinal::PositionDataFrame(data.frame(x = 1:2, y = 1)), centroided = FALSE)
        actual <- suppressMessages(estimateFwhm(profile))
        expect_equal(getFwhm(actual, f$predictionMasses), getFwhm(model, f$predictionMasses), tolerance = 1e-10)
        # This is a generated PROFILE pair for adapter testing, not the supplied
        # reduced centroided acquisition used by extraction/map benchmarks.
        path <- write_msi_pair(profile, tempfile(fileext = ".imzML"))
        backed <- readMsi(path, verbose = FALSE)
        expect_equal(getFwhm(suppressMessages(estimateFwhm(backed)), f$predictionMasses),
                     getFwhm(model, f$predictionMasses), tolerance = 1e-10)
    }
    table <- data.frame(mz = f$widths$peaks, fwhm = f$widths$fwhmValues, snr = 10)
    expect_equal(getFwhm(estimateFwhm(table), f$predictionMasses), f$predictions, tolerance = 1e-10)
    smooth <- estimateFwhm(table, smoothingMethod = "superSmoother")
    expect_true(all(is.finite(getFwhm(smooth, f$predictionMasses))))
    expect_error(estimateFwhm(make_msi_arrays()), "profile spectra")
    expect_error(estimateFwhm(make_msi_arrays(centroided = NA)), "profile spectra")
    multi <- make_msi_arrays(rep(list(f$profile$mz), 2), rep(list(f$profile$intensity), 2),
        data.frame(x = 1:2, y = 1), runs = c("a", "b"), centroided = FALSE)
    expect_error(estimateFwhm(multi), "one run")
    expect_error(estimateFwhm(Cardinal::normalize(make_msi_arrays())), "Cardinal::process")
    table$fwhm[1] <- -1
    expect_error(estimateFwhm(table), "Precomputed fwhm")
    for (bad in c(0, NA_real_, Inf, 1.5)) {
        expect_error(estimateFwhm(f$profile, peakSampling = bad), "peakSampling")
    }
    expect_error(suppressMessages(estimateFwhm(list(mz = 1:20, intensity = rep(1, 20)))), "five distinct")
})

test_that("native TIC and RMS use Cardinal scaling and can request division-only results", {
    f <- phase4_fixture()
    for (arrays in c(TRUE, FALSE)) {
        x <- phase4_input(arrays)
        for (method in c("TIC", "RMS")) {
            actual <- normalizeIntensity(x, method)
            direct <- Cardinal::process(Cardinal::normalize(x, tolower(method)),
                verbose = FALSE, BPPARAM = BiocParallel::SerialParam())
            expect_equal(phase4_values(actual), phase4_values(direct), tolerance = 1e-10,
                         ignore_attr = TRUE)
            expect_true(methods::is(Cardinal::intensity(actual), "matter"))
            expect_length(ProtGenerics::processingData(actual), 0)
            expect_true(length(attr(ProtGenerics::processingData(actual), "processed")) >= 1)
            expect_equal(.msiCoordinates(actual), .msiCoordinates(x))
            expect_equal(.msiRun(actual), .msiRun(x))
            again <- normalizeIntensity(actual, method)
            expect_equal(phase4_values(again), phase4_values(actual), tolerance = 1e-10, ignore_attr = TRUE)
            if (arrays) {
                division <- normalizeIntensity(x, method, scale = 1)
                expect_equal(phase4_values(division), f$normalizations[[method]],
                             tolerance = 1e-10, ignore_attr = TRUE)
            }
        }
    }
})

test_that("reference normalization follows Cardinal zero-reference policy within explicit tolerance", {
    x <- phase4_input(FALSE)
    actual <- normalizeIntensity(x, "IS", mz = 200, tol = 0)
    direct <- Cardinal::process(Cardinal::normalize(x, "reference", ref = 200),
        verbose = FALSE, BPPARAM = BiocParallel::SerialParam())
    expect_equal(phase4_values(actual), phase4_values(direct), tolerance = 1e-10, ignore_attr = TRUE)
    expect_equal(phase4_values(actual)[, 2], numeric(5), ignore_attr = TRUE)
    arrays <- normalizeIntensity(phase4_input(), "IS", mz = 200, tol = 0)
    expect_equal(phase4_values(arrays)[[2]], numeric(4))
    expect_true(is.na(Cardinal::pixelData(arrays)$moleculaRNormalizationFactor[2]))
    # An outside standard cannot silently use an unrelated nearest mass.
    absent <- normalizeIntensity(phase4_input(), "IS", mz = 500, tol = 0.1)
    expect_true(all(unlist(phase4_values(absent)) == 0))
    expect_error(normalizeIntensity(x, "IS", mz = 200), "tol")
})

test_that("median and aligned fold-change workers preserve documented observed/zero semantics", {
    f <- phase4_fixture()
    expect_equal(phase4_values(normalizeIntensity(phase4_input(), "median")),
                 f$normalizations$median, tolerance = 1e-10, ignore_attr = TRUE)
    for (method in c("mean", "median", "sum")) {
        normalized <- normalizeIntensity(phase4_input(FALSE), "medianFoldChange", meanMethod = method,
                                         featureChunkSize = 2)
        values <- phase4_values(normalized)
        for (i in 1:3) expect_equal(unname(values[f$signal[, i] > 0, i]), f$fold[[method]][[i]], tolerance = 1e-10)
    }
    stored <- normalizeIntensity(phase4_input(FALSE), "medianFoldChange", zeroAsMissing = FALSE)
    expect_equal(lapply(1:3, function(i) unname(phase4_values(stored)[, i])), f$fullFold, tolerance = 1e-10)
    expect_error(normalizeIntensity(phase4_input(), "medianFoldChange"), "aligned")
    invalid <- make_msi_experiment(c(100, 200, 300), matrix(c(0, 0, 1), ncol = 1), data.frame(x = 1, y = 1))
    expect_error(normalizeIntensity(invalid, "median"), "factor")
})

test_that("normalization validates queues and output, keeps zero spectra, and preserves source files", {
    x <- make_msi_arrays(list(100, numeric()), list(0, numeric()), data.frame(x = 1:2, y = 1))
    for (method in c("TIC", "RMS", "median", "IS")) {
        actual <- normalizeIntensity(x, method, mz = 100, tol = 0)
        expect_equal(phase4_values(actual), list(0, numeric()))
    }
    expect_error(normalizeIntensity(Cardinal::normalize(phase4_input()), "TIC"), "Cardinal::process")
    expect_error(normalizeIntensity(phase4_input(), "bad"))
    expect_error(normalizeIntensity(phase4_input(), "TIC", scale = 0), "scale")
    for (bad in c(NA_real_, NaN, Inf, -1)) {
        invalid <- make_msi_arrays(list(100), list(bad), data.frame(x = 1, y = 1))
        expect_error(normalizeIntensity(invalid, "TIC"), "intensities")
    }
    path <- write_msi_pair(phase4_input(), tempfile(fileext = ".imzML"))
    inputFiles <- c(path, paste0(tools::file_path_sans_ext(path), ".ibd"))
    before <- tools::md5sum(inputFiles)
    backed <- readMsi(path, verbose = FALSE)
    result <- normalizeIntensity(backed, "TIC")
    expect_identical(tools::md5sum(inputFiles), before)
    expect_error(normalizeIntensity(backed, "TIC", outfile = path), "already exists")
    expect_true(all(file.exists(.msiBackingState(result)$path)))
    cache <- new.env(parent = emptyenv())
    original <- .msiExtractIons(backed, 100, 0.1, cache = cache)
    expect_false(identical(as.numeric(.msiExtractIons(result, 100, 0.1, cache = cache)), as.numeric(original)))
})

test_that("frequency filtering counts positive pixels inclusively and preserves metadata", {
    x <- phase4_input(FALSE)
    actual <- filterPeaks(x, minFreq = 1, featureChunkSize = 2)
    expect_equal(Cardinal::mz(actual), c(99.9, 100, 100.1))
    expect_equal(as.numeric(Cardinal::featureData(actual)$count), rep(3, 3))
    expect_equal(as.numeric(Cardinal::featureData(actual)$freq), rep(1, 3))
    equal <- filterPeaks(x, minFreq = 2/3)
    expect_equal(Cardinal::mz(equal), Cardinal::mz(x))
    expect_equal(as.numeric(Cardinal::featureData(equal)$count), c(3, 3, 3, 2, 2))
    expect_equal(.msiCoordinates(actual), .msiCoordinates(x))
    expect_equal(.msiRun(actual), .msiRun(x))
    expect_equal(Cardinal::experimentData(actual), Cardinal::experimentData(x))
    expect_length(ProtGenerics::processingData(actual), 0)
    duplicate <- make_msi_experiment(c(100, 100), matrix(c(1, 0, 0, 1, 0, 0), nrow = 2),
                                     data.frame(x = 1:3, y = 1))
    expect_equal(nrow(filterPeaks(duplicate, 1/3)), 2)
    expect_equal(nrow(filterPeaks(duplicate, 0.5)), 0)
    empty <- make_msi_experiment(c(100, 200), matrix(0, 2, 3), data.frame(x = 1:3, y = 1))
    removed <- filterPeaks(empty, 0.1)
    expect_equal(nrow(removed), 0)
    expect_equal(.msiExtractIon(removed, 100, 0.1), numeric(3))
    expect_equal(nrow(filterPeaks(empty, 0)), 2)
    expect_error(filterPeaks(phase4_input(), 0.5), "aligned")
    expect_error(filterPeaks(Cardinal::normalize(x), 0.5), "Cardinal::process")
    for (bad in c(NA_real_, Inf, -1, 1.1)) expect_error(filterPeaks(x, bad), "minFreq")
    invalid <- make_msi_experiment(c(100, 200), matrix(c(NA, 1, 0, 1, 0, 1), 2, 3),
                                   data.frame(x = 1:3, y = 1))
    expect_error(filterPeaks(invalid, 0.5), "intensities")
    expect_error(filterPeaks(x, 0.5, featureChunkSize = 0), "featureChunkSize")
})

test_that("native preprocessing reaches seeded maps without MALDIquant algorithms", {
    f <- legacy_fixture()
    profile <- phase4_fixture()
    local_mocked_bindings(detectPeaks = function(...) stop("legacy picker called"),
        createMassSpectrum = function(...) stop("legacy spectrum created"), .package = "MALDIquant")
    model <- suppressMessages(estimateFwhm(profile$profile))
    x <- make_msi_experiment(f$mz, f$intensity, f$coordinates)
    normalized <- normalizeIntensity(x, "TIC")
    division <- normalizeIntensity(x, "TIC", scale = 1)
    filtered <- filterPeaks(normalized, 0.01)
    batch <- batchAnalyteSearch(filtered, model, m = 100, verbose = FALSE)
    single <- searchAnalyte(100, getFwhm(model, 100), division)
    expect_s3_class(batch, "analytePointPattern")
    expect_equal(batch$marks$intensity, 5 * single$marks$intensity, tolerance = 1e-10)
    first <- suppressWarnings(probMap(batch, bw = 1, seed = 1, verbose = FALSE))
    second <- suppressWarnings(probMap(single, bw = 1, seed = 1, verbose = FALSE))
    expect_equal(first$rhoMoi$v, second$rhoMoi$v, tolerance = 1e-10)
    expect_equal(first$rhoCsr$v, second$rhoCsr$v, tolerance = 1e-10)
    expect_equal(first$hotspotIm$v, second$hotspotIm$v, tolerance = 1e-10)
    expect_equal(first$coldspotIm$v, second$coldspotIm$v, tolerance = 1e-10)
})

test_that("native preprocessing workers preserve library environment and backend ownership", {
    backend <- BiocParallel::SnowParam(workers = 2L, progressbar = FALSE)
    previous <- Sys.getenv("R_LIBS", unset = NA_character_)
    on.exit({
        BiocParallel::bpstop(backend)
        if (is.na(previous)) Sys.unsetenv("R_LIBS") else Sys.setenv(R_LIBS = previous)
    }, add = TRUE)
    Sys.unsetenv("R_LIBS")
    actual <- normalizeIntensity(phase4_input(), "TIC", BPPARAM = backend)
    expect_true(is.na(Sys.getenv("R_LIBS", unset = NA_character_)))
    expect_false(BiocParallel::bpisup(backend))
    expected <- normalizeIntensity(phase4_input(), "TIC")
    expect_equal(phase4_values(actual), phase4_values(expected), tolerance = 1e-10, ignore_attr = TRUE)
    .msiStartBackend(backend)
    actual <- normalizeIntensity(phase4_input(), "RMS", BPPARAM = backend)
    expect_true(BiocParallel::bpisup(backend))
    expect_equal(phase4_values(actual), phase4_values(normalizeIntensity(phase4_input(), "RMS")),
                 tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("Cardinal byte budgets remain distinct from record counts", {
    previous <- Cardinal::getCardinalChunksize()
    on.exit(Cardinal::setCardinalChunksize(previous), add = TRUE)
    Cardinal::setCardinalChunksize(c(B = 48))
    expect_equal(lengths(.msiFeatureBlocks(phase4_input(FALSE))), c(2L, 2L, 1L))
    expect_equal(.msiBatchChunkSize(3, NULL, make_msi_arrays()), 1)
    expect_equal(.msiBatchChunkSize(3, 2, make_msi_arrays()), 2)
    expect_equal(.msiExtractIons(make_msi_arrays(), 100, 0.1)[, 1], c(10.375, 3, 0))
    Cardinal::setCardinalChunksize(2)
    expect_equal(.msiBatchChunkSize(3, NULL, make_msi_arrays()), 2)
    expect_equal(lengths(.msiFeatureBlocks(phase4_input(FALSE))), 5L)
})
