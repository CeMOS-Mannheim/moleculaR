# Batch references are independent repeated single searches, including both
# representation-specific means and immutable pre-migration spatial fixtures.
batch_model <- function(width = 0.1) fwhm(function(mz) rep(width, length(mz)))
batch_window <- function() spatstat.geom::owin(c(0, 4), c(0, 3))

batch_fields <- function(spp) {
    if (!spp$n) return(list())
    lapply(seq_len(nrow(spp$metaData)), function(i) {
        points <- which(spp$marks$idx == spp$metaData$idx[i])
        list(x = spp$x[points], y = spp$y[points], intensity = spp$marks$intensity[points],
             metadata = spp$metaData[i, setdiff(names(spp$metaData), "idx"), drop = FALSE])
    })
}

test_that("chunked sparse extraction equals independent single queries", {
    masses <- c(101, 100, 500, 99.9, 100)
    widths <- c(0.1, 0.1, 0.1, 0.01, 0.1)
    for (x in list(make_msi_arrays(), make_msi_experiment())) {
        for (method in c("Gaussian", "sum", "mean", "max")) {
            expected <- vapply(seq_along(masses), function(i) {
                .msiExtractIon(x, masses[i], widths[i], method)
            }, numeric(3))
            for (size in c(1, 2, 20)) {
                actual <- .msiExtractIons(x, masses, widths, method, chunkSize = size)
                expect_s4_class(actual, "dgCMatrix")
                expect_equal(as.matrix(actual), expected, tolerance = 1e-10)
            }
        }
        expect_equal(dim(.msiExtractIons(x, numeric(), numeric())), c(3L, 0L))
        expect_equal(as.matrix(.msiExtractIons(x, 500, 0.1)), matrix(0, 3, 1))
    }
})

test_that("batch boundaries, empty spectra and numeric validation preserve single semantics", {
    width <- 2 * sqrt(2 * log(2))
    x <- make_msi_arrays(list(c(104, 103, 100, 97, 100, 96), numeric()),
        list(c(1e12, 2, 3, 4, 5, 1e12), numeric()), data.frame(x = 1:2, y = 1))
    for (method in c("Gaussian", "sum", "mean", "max")) {
        expected <- .msiExtractIon(x, 100, width, method)
        expect_equal(as.numeric(.msiExtractIons(x, 100, width, method)), expected)
    }
    for (bad in list(0, -1, Inf, NaN, NA_real_, "100")) {
        expect_error(.msiExtractIons(x, bad, 0.1))
        expect_error(.msiExtractIons(x, 100, bad))
    }
    expect_error(.msiExtractIons(x, c(100, 101), 0.1), "equal length")
    expect_error(.msiExtractIons(x, 100, 0.1, "bad"), "wMethod")
    for (bad in list(0, -1, NA_real_, Inf, 1.5, "2")) {
        expect_error(.msiExtractIons(x, 100, 0.1, chunkSize = bad), "chunkSize")
    }
    expect_error(.msiExtractIons(x, 100, 0.1, BPPARAM = 1), "BPPARAM")
    expect_error(.msiExtractIons(x, 100, 0.1, cache = list()), "cache")
    expect_error(.msiExtractIons(x, 100, 0.1, verbose = NA), "verbose")
    for (bad in c(NA_real_, NaN, Inf, -1)) {
        invalid <- make_msi_arrays(list(100), list(bad), data.frame(x = 1, y = 1))
        expect_error(.msiExtractIons(invalid, 100, 0.1), "intensities")
    }
    expect_error(.msiAggregateIons(c(1e308, 1e308), c(100, 100),
        .msiBatchQueries(100, 0.1, "sum"), "sum"), "not finite")
    expect_equal(.msiAggregateIons(1e-300, 100, .msiBatchQueries(100, 0.1, "sum"), "sum")$x, 1e-300)
})

test_that("public batch output retains target order, duplicates, confirmations and acquired domain", {
    x <- make_msi_arrays()
    masses <- c(101, NA, 100, 500, 100, 99.9)
    for (method in c("Gaussian", "sum", "mean", "max")) {
        actual <- batchAnalyteSearch(x, batch_model(), batch_window(), masses,
            wMethod = method, verifiedMasses = c(100, 101), verbose = FALSE, chunkSize = 2)
        single <- lapply(masses[!is.na(masses)], function(mass) {
            searchAnalyte(mass, 0.1, x, wMethod = method, spwin = batch_window(),
                verifiedMasses = c(100, 101))
        })
        expected <- superimposeAnalytes(single, spWin = batch_window(), check = FALSE)
        expect_equal(batch_fields(actual), batch_fields(expected), ignore_attr = TRUE)
        expect_true(all(actual$marks$idx %in% actual$metaData$idx))
        expect_equal(spatstat.geom::area.owin(actual$window), 12)
    }
    confirmed <- batchAnalyteSearch(x, batch_model(), batch_window(), c(100, 101),
        verifiedMasses = 100, confirmedOnly = TRUE, verbose = FALSE)
    expect_equal(confirmed$metaData$mzVals, 100)
    for (empty in list(numeric(), NA_real_, 500)) {
        expect_equal(batchAnalyteSearch(x, batch_model(), batch_window(), empty, verbose = FALSE)$n, 0)
    }
    # The query centre can be outside the mass range while its interval overlaps.
    edge <- batchAnalyteSearch(x, batch_model(), batch_window(), 99.85, verbose = FALSE)
    expect_equal(batch_fields(edge), batch_fields(searchAnalyte(99.85, 0.1, x,
        spwin = batch_window())), ignore_attr = TRUE)
    expect_equal(batchAnalyteSearch(x, batch_model(), batch_window(), 100,
        confirmedOnly = TRUE, verbose = FALSE)$n, 0)
    for (bad in list("100", 0, Inf, NaN)) {
        expect_error(batchAnalyteSearch(x, batch_model(), batch_window(), bad, verbose = FALSE), "m must")
    }
    expect_error(batchAnalyteSearch(x, batch_model(NA_real_), batch_window(), 100, verbose = FALSE), "fwhm")
    expect_error(batchAnalyteSearch(x, batch_model(), batch_window(), 100,
        confirmedOnly = NA, verbose = FALSE), "confirmedOnly")
    expect_error(batchAnalyteSearch(x, batch_model(), batch_window(), numeric(),
        verifiedMasses = "100", verbose = FALSE), "m/z")
    expect_error(.msiBatchBackend(0, NULL), "numCores")
    expect_error(.msiBatchBackend(2, BiocParallel::SerialParam()), "either")
    expect_s4_class(.msiBatchBackend(2, NULL), "SnowParam")
    expect_error(batchAnalyteSearch(Cardinal::normalize(x), batch_model(), batch_window(),
        100, verbose = FALSE), "Cardinal::process")
})

batch_lipids <- function() {
    data.frame("Lipid ID" = c("lipid-a", "lipid-b"), "Formula (pH7.3)" = c("A", "B"),
        "Abbreviation*" = c("PC(30:1)", "PC(32:2)"), numDoubleBond = 1:2,
        lipidGroup = "PC(x:x)", chainLength = c(30L, 32L),
        "Exact m/z of [M-H]-" = c(100, 100), "Exact m/z of [M+H]+" = c(101, NA),
        "Exact m/z of [M+Na]+" = c(99.9, 500), "Exact m/z of [M+K]+" = c(100.1, 100),
        check.names = FALSE)
}

test_that("lipid queries use one engine while preserving adduct and database metadata", {
    x <- make_msi_arrays()
    database <- batch_lipids()
    queries <- .msiLipidQueries(database, c("M+K", "M-H", "M+H", "M+Na"))
    expect_equal(queries$mz, c(100, 101, 99.9, 100.1, 100, NA, 500, 100))
    for (method in c("Gaussian", "sum", "mean", "max")) {
        actual <- batchLipidSearch(x, batch_model(), batch_window(), database,
            wMethod = method, verifiedMasses = 100, verbose = FALSE)
        singles <- lapply(which(!is.na(queries$mz)), function(i) {
            searchAnalyte(queries$mz[i], 0.1, x, wMethod = method,
                spwin = batch_window(), verifiedMasses = 100, metaData = queries$metadata[[i]])
        })
        expected <- superimposeAnalytes(singles, spWin = batch_window(), check = FALSE)
        expect_equal(batch_fields(actual), batch_fields(expected), ignore_attr = TRUE)
        expect_true(all(c("mode", "adduct", "lipidID", "sumformula", "abbrev",
            "numDoubleBonds", "lipidClass", "chainLength") %in% names(actual$metaData)))
        expect_equal(actual$metaData$lipidID, c(rep("lipid-a", 4), rep("lipid-b", 2)))
        expect_equal(actual$metaData$mode, c("negative", "positive", "positive", "positive",
                                            "negative", "positive"))
        expect_equal(actual$metaData$adduct, c("M-H", "M+H", "M+Na", "M+K", "M-H", "M+K"))
        expect_equal(actual$metaData$sumformula, c(rep("A", 4), rep("B", 2)))
        expect_equal(actual$metaData$numDoubleBonds, c(rep(1L, 4), rep(2L, 2)))
        expect_equal(actual$metaData$chainLength, c(rep(30L, 4), rep(32L, 2)))
    }
    for (adduct in c("M-H", "M+H", "M+Na", "M+K", "M+k")) {
        actual <- batchLipidSearch(x, batch_model(), batch_window(), database,
            adduct = adduct, verbose = FALSE)
        expect_true(all(actual$metaData$adduct == sub("M+k", "M+K", adduct, fixed = TRUE)))
    }
    # 0.12 is outside one FWHM but inside three sigma; the old prefilter dropped it.
    confirmed <- batchLipidSearch(x, batch_model(), batch_window(), database,
        adduct = "M-H", verifiedMasses = 100.12, confirmedOnly = TRUE, verbose = FALSE)
    expect_equal(confirmed$metaData$lipidID, database[["Lipid ID"]])
    expect_true(all(confirmed$metaData$mzConfirmed))
    expect_equal(batchLipidSearch(x, batch_model(), batch_window(), database,
        verifiedMasses = numeric(), confirmedOnly = TRUE, verbose = FALSE)$n, 0)
    expect_equal(batchLipidSearch(x, batch_model(), batch_window(), database[FALSE, ], verbose = FALSE)$n, 0)
    expect_equal(batchLipidSearch(x, batch_model(), batch_window(), database,
        adduct = character(), verbose = FALSE)$n, 0)
    expect_error(.msiLipidQueries(database, "bad"), "adduct")
    expect_error(.msiLipidQueries(database[, names(database) != "lipidGroup"], "M-H"), "sldb")
    database[["Exact m/z of [M-H]-"]] <- "100"
    expect_error(.msiLipidQueries(database, "M-H"), "numeric")
})

test_that("file-backed chunks and batch-derived probability maps match single searches", {
    for (arrays in c(TRUE, FALSE)) {
        path <- write_msi_pair(fixture_cardinal(arrays), tempfile(fileext = ".imzML"))
        x <- readMsi(path, as = if (arrays) "MSImagingArrays" else "MSImagingExperiment", verbose = FALSE)
        expect_true(methods::is(Cardinal::intensity(x), "matter"))
        actual <- batchAnalyteSearch(x, batch_model(), m = c(100, 500), verbose = FALSE, chunkSize = 40)
        expected <- searchAnalyte(100, 0.1, x)
        expect_equal(batch_fields(actual), batch_fields(expected), ignore_attr = TRUE)
        map <- suppressWarnings(probMap(actual, bw = 1, seed = 1, verbose = FALSE))
        reference <- legacy_fixture()$mpm
        expect_equal(map$rhoMoi$v, reference$rhoMoi, tolerance = 1e-10)
        expect_equal(map$rhoCsr$v, reference$rhoCsr, tolerance = 1e-10)
        expect_equal(map$hotspotIm$v, reference$hotspot, tolerance = 1e-10)
        expect_equal(map$coldspotIm$v, reference$coldspot, tolerance = 1e-10)
    }
})

test_that("optional cache reuses only equivalent spectra, queries and processing states", {
    x <- make_msi_arrays()
    cache <- new.env(parent = emptyenv())
    first <- .msiExtractIons(x, c(100, 101), c(0.1, 0.1), cache = cache)
    expect_identical(.msiExtractIons(x, c(100, 101), c(0.1, 0.1), cache = cache), first)
    .msiExtractIons(x, 100, 0.2, "sum", cache = cache)
    expect_identical(cache$.moleculaR_ions$key$fwhm, 0.2)
    expect_identical(cache$.moleculaR_ions$key$method, "sum")
    normalized <- Cardinal::process(Cardinal::normalize(x, method = "tic"),
        verbose = FALSE, BPPARAM = BiocParallel::SerialParam())
    changed <- .msiExtractIons(normalized, 100, 0.1, cache = cache)
    expect_false(identical(as.numeric(changed), as.numeric(first[, 1])))
    expect_identical(cache$.moleculaR_ions$key$object, normalized)
    roi <- .msiSelectPixels(x, c(2, 1))
    expect_equal(nrow(.msiExtractIons(roi, 100, 0.1, cache = cache)), 2)
    path <- write_msi_pair(make_msi_arrays(), tempfile(fileext = ".imzML"))
    backed <- readMsi(path, verbose = FALSE)
    .msiExtractIons(backed, 100, 0.1, cache = cache)
    previous <- cache$.moleculaR_ions$key$files
    expect_true(nrow(previous) >= 1)
    Sys.setFileTime(previous$path[1], as.POSIXct(previous$mtime[1], origin = "1970-01-01") + 10)
    .msiExtractIons(backed, 100, 0.1, cache = cache)
    expect_false(identical(cache$.moleculaR_ions$key$files, previous))
    expect_length(ls(cache, all.names = TRUE), 1)
})

test_that("a batch reads each spectrum once and cache hits perform no spectral reads", {
    x <- make_msi_arrays()
    calls <- 0L
    aggregate <- .msiAggregateIons
    local_mocked_bindings(.msiAggregateIons = function(...) {
        calls <<- calls + 1L
        aggregate(...)
    })
    cache <- new.env(parent = emptyenv())
    masses <- c(101, 100, 99.9, 500, 100)
    widths <- rep(0.1, length(masses))
    .msiExtractIons(x, masses, widths, chunkSize = 2, cache = cache)
    expect_equal(calls, 3L)
    .msiExtractIons(x, masses, widths, chunkSize = 1, cache = cache)
    expect_equal(calls, 3L)
    .msiExtractIons(x, masses, widths * 2, cache = cache)
    expect_equal(calls, 6L)
    .msiExtractIons(x, masses, widths * 2, "sum", cache = cache)
    expect_equal(calls, 9L)
    changed <- make_msi_arrays(intensity = list(c(2, 20, 4), c(3, 9), numeric()))
    .msiExtractIons(changed, masses, widths * 2, "sum", cache = cache)
    expect_equal(calls, 12L)
    .msiExtractIons(make_msi_experiment(), masses, widths * 2, "sum", cache = cache)
    expect_equal(calls, 15L)
})

test_that("Snow workers agree with serial extraction and respect backend ownership", {
    backend <- BiocParallel::SnowParam(workers = 2L, progressbar = FALSE)
    on.exit(BiocParallel::bpstop(backend), add = TRUE)
    expect_false(BiocParallel::bpisup(backend))
    x <- make_msi_arrays()
    expected <- .msiExtractIons(x, c(100, 101), c(0.1, 0.1))
    expect_equal(.msiExtractIons(x, c(100, 101), c(0.1, 0.1), BPPARAM = backend), expected)
    expect_false(BiocParallel::bpisup(backend))
    BiocParallel::bpstart(backend)
    for (arrays in c(TRUE, FALSE)) {
        native <- if (arrays) make_msi_arrays() else make_msi_experiment()
        path <- write_msi_pair(native, tempfile(fileext = ".imzML"))
        backed <- readMsi(path, as = if (arrays) "MSImagingArrays" else "MSImagingExperiment", verbose = FALSE)
        for (method in c("Gaussian", "sum", "mean", "max")) {
            expected <- .msiExtractIons(backed, c(100, 101), c(0.1, 0.1), method)
            expect_equal(.msiExtractIons(backed, c(100, 101), c(0.1, 0.1), method,
                chunkSize = 3, BPPARAM = backend), expected, tolerance = 1e-10)
        }
    }
    expect_true(BiocParallel::bpisup(backend))
    expected <- batchAnalyteSearch(x, batch_model(), batch_window(), c(100, 101), verbose = FALSE)
    parallel <- batchAnalyteSearch(x, batch_model(), batch_window(), c(100, 101),
        numCores = 2, verbose = FALSE)
    expect_equal(batch_fields(parallel), batch_fields(expected), ignore_attr = TRUE)
})

test_that("batch compatibility path and default chunk configuration remain usable", {
    f <- legacy_fixture()
    sparse <- sparseIntensityMatrix(Matrix::Matrix(t(f$intensity), sparse = TRUE), f$mz, f$coordinates)
    # The retained legacy window diagnostic warns on the fixture's hole.
    actual <- suppressWarnings(batchAnalyteSearch(sparse, batch_model(),
        m = c(100, 100), verbose = FALSE))
    expect_equal(actual$n, 2 * length(f$patterns$Gaussian$x))
    expect_error(batchAnalyteSearch(sparse, batch_model(), m = 100,
        cache = new.env(), verbose = FALSE), "Cardinal input")
    expect_error(batchLipidSearch(sparse, batch_model(), sldb = batch_lipids(),
        BPPARAM = BiocParallel::SerialParam()), "Cardinal input")
    previous <- Cardinal::getCardinalChunksize()
    on.exit(Cardinal::setCardinalChunksize(previous), add = TRUE)
    Cardinal::setCardinalChunksize(2)
    expect_equal(.msiBatchChunkSize(3, NULL), 2)
    expect_equal(.msiBatchChunkSize(3, 1), 1)
})

test_that("cache hits refresh lipid labels and confirmation without re-reading spectra", {
    x <- make_msi_arrays()
    database <- batch_lipids()
    cache <- new.env(parent = emptyenv())
    calls <- 0L
    aggregate <- .msiAggregateIons
    local_mocked_bindings(.msiAggregateIons = function(...) {
        calls <<- calls + 1L
        aggregate(...)
    })
    first <- batchLipidSearch(x, batch_model(), batch_window(), database,
        verifiedMasses = 100, verbose = FALSE, cache = cache)
    expect_equal(calls, 3L)
    database[["Lipid ID"]] <- c("new-a", "new-b")
    second <- batchLipidSearch(x, batch_model(), batch_window(), database,
        verifiedMasses = 101, verbose = FALSE, cache = cache)
    expect_equal(calls, 3L)
    expect_equal(second$metaData$lipidID, c(rep("new-a", 4), rep("new-b", 2)))
    expect_equal(second$metaData$mzConfirmed, second$metaData$mzVals == 101)
    expect_equal(first$marks$intensity, second$marks$intensity)
})
