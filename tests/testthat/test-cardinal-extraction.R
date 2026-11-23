test_that("all aggregations follow hand-calculated spectra in pixel order", {
    arrays <- make_msi_arrays()
    experiment <- make_msi_experiment()
    expected <- list(Gaussian = c(10.375, 3, 0), sum = c(16, 3, 0),
                     max = c(10, 3, 0))
    for (method in names(expected)) {
        expect_equal(.msiExtractIon(arrays, 100, 0.1, method), expected[[method]], tolerance = 1e-10)
        expect_equal(.msiExtractIon(experiment, 100, 0.1, method), expected[[method]], tolerance = 1e-10)
    }
    expect_equal(.msiExtractIon(arrays, 100, 0.1, "mean"), c(16 / 3, 3, 0))
    expect_equal(.msiExtractIon(experiment, 100, 0.1, "mean"), c(16 / 3, 1, 0))
    for (x in list(arrays, experiment)) {
        expect_equal(.msiExtractIon(x, 500, 0.1), numeric(3))
        expect_equal(.msiExtractIon(x, 99.9, 1e-8, "sum"), c(2, 0, 0))
    }
})

test_that("inclusive boundaries, duplicate and unsorted peaks are explicit", {
    width <- 2 * sqrt(2 * log(2)) # sigma = 1; exact integer boundaries
    x <- make_msi_arrays(mz = list(c(104, 103, 100, 97, 100, 96)),
                        intensity = list(c(100, 2, 3, 4, 5, 200)),
                        coordinates = data.frame(x = 1, y = 1))
    expect_equal(.msiExtractIon(x, 100, width, "sum"), 14)
    expect_equal(.msiExtractIon(x, 100, width), 8 + 6 * exp(-4.5))
    expect_equal(.msiExtractIon(x, 100, width, "mean"), 3.5)
    expect_equal(.msiExtractIon(x, 100, width, "max"), 5)
    for (method in c("Gaussian", "sum", "mean", "max")) {
        empty <- make_msi_arrays(list(numeric()), list(numeric()), data.frame(x = 1, y = 1))
        expect_equal(.msiExtractIon(empty, 100, 0.1, method), 0)
    }
})

test_that("invalid targets and retrieved spectral values fail clearly", {
    x <- make_msi_arrays()
    for (bad in list(numeric(), c(1, 2), NA_real_, NaN, Inf, 0, -1, "100")) {
        expect_error(.msiExtractIon(x, bad, 0.1), "m must")
        expect_error(.msiExtractIon(x, 100, bad), "fwhm must")
    }
    for (bad in list("bad", NA_character_, c("sum", "mean")))
        expect_error(.msiExtractIon(x, 100, 0.1, bad), "wMethod")
    for (bad in c(NA_real_, NaN, Inf, -1)) {
        expect_error(.msiAggregateIon(c(1, bad), c(100, 101), 100, 0.1, "sum", 0.2), "intensities")
        expect_error(.msiAggregateIon(c(1, 2), c(100, bad), 100, 0.1, "sum", 0.2), "m/z")
    }
    expect_error(.msiAggregateIon(c(1, 2), 100, 100, 0.1, "sum", 0.2), "match")
    expect_error(.msiAggregateIon(c(1e308, 1e308), c(100, 100), 100, 0.1, "sum", 0.2), "not finite")
    expect_equal(.msiAggregateIon(1e-300, 100, 100, 0.1, "sum", 0.2), 1e-300)
    expect_error(searchAnalyte(100, 0.1, list()), "spData")
})

test_that("Cardinal metadata, coordinates, run selection and processing stay valid", {
    x <- make_msi_arrays()
    expect_true(.msiIsObject(x))
    expect_false(.msiIsObject(list()))
    expect_equal(.msiCoordinates(x), data.frame(x = c(3, 1, 2), y = c(1, 2, 1)))
    expect_equal(.msiNumberOfPixels(x), 3L)
    expect_equal(.msiRun(x), rep("sample", 3))
    expect_error(.msiValidate(list()), "Cardinal")
    for (status in c(FALSE, NA))
        expect_error(.msiExtractIon(make_msi_arrays(centroided = status), 100, 0.1), "centroided")
    queued <- Cardinal::normalize(x, method = "tic")
    expect_error(.msiExtractIon(queued, 100, 0.1), "Cardinal::process")
    processed <- Cardinal::process(queued, verbose = FALSE, BPPARAM = BiocParallel::SerialParam())
    expect_equal(.msiExtractIon(processed, 100, 0.1, "sum"), c(3, 0.75, 0))
    expect_length(ProtGenerics::processingData(queued), 1)
    expect_length(ProtGenerics::processingData(processed), 0)
    multi <- make_msi_arrays(runs = c("a", "b", "a"))
    expect_silent(.msiValidate(multi))
    expect_error(.msiExtractIon(multi, 100, 0.1), "single run")
    subset <- .msiSelectPixels(multi, .msiRun(multi) == "b")
    expect_equal(.msiExtractIon(subset, 100, 0.1, "sum"), 3)
    for (object in list(x, make_msi_experiment())) {
        expect_equal(.msiCoordinates(.msiSelectPixels(object, c(3, 1))),
                     .msiCoordinates(object)[c(3, 1), , drop = FALSE], ignore_attr = TRUE)
        for (bad in list(integer(), 0, 4, 1.5, c(1, 1), NA_integer_, c(TRUE, FALSE)))
            expect_error(.msiSelectPixels(object, bad), "pixel")
    }
    expect_error(.msiValidate(make_msi_arrays(coordinates = data.frame(x = c(1, 1, 2), y = c(1, 1, 2)))), "Duplicate")
    expect_error(.msiValidate(make_msi_arrays(coordinates = data.frame(x = c(1, Inf, 2), y = c(1, 1, 2)))), "finite")
    expect_error(.msiValidate(make_msi_arrays(coordinates = data.frame(x = 1:3, y = 1:3, z = 1:3))), "z plane")
})

test_that("annotation matching uses absolute tolerance with deterministic ties", {
    expect_equal(.msiNearestMass(100, c(100.1, 99.9), 0.2), 1L)
    expect_equal(.msiNearestMass(100, c(110, 100.01, 100), 0.1), 3L)
    expect_true(is.na(.msiNearestMass(100, numeric(), 0.1)))
    expect_true(is.na(.msiNearestMass(100, NA, 0.1)))
    expect_true(is.na(.msiNearestMass(100, c(99, 101), 0.1)))
    for (bad in list("100", c(NA, 100), Inf))
        expect_error(.msiNearestMass(100, bad, 0.1), "m/z")
    win <- spatstat.geom::owin(c(0, 4), c(0, 3))
    x <- make_msi_arrays()
    yes <- searchAnalyte(100, 0.1, x, spwin = win, verifiedMasses = c(110, 100),
                         confirmedOnly = TRUE, metaData = list(label = "known"))
    expect_true(yes$metaData$mzConfirmed)
    expect_equal(yes$metaData$label, "known")
    expect_true(all(yes$marks$idx %in% yes$metaData$idx))
    expect_equal(yes$x, c(3, 1))
    expect_equal(yes$y, c(1, 2))
    expect_equal(searchAnalyte(100, 0.1, x, spwin = win, confirmedOnly = TRUE)$n, 0)
    expect_false(searchAnalyte(100, 0.1, x, spwin = win, verifiedMasses = numeric())$metaData$mzConfirmed)
    expect_error(searchAnalyte(100, 0.1, x, spwin = win, confirmedOnly = NA), "confirmedOnly")
    expect_error(searchAnalyte(100, 0.1, x, spwin = win, metaData = list(idx = 1)), "metaData")
    expect_error(searchAnalyte(100, 0.1, x, spwin = "bad"), "spwin")
    expect_error(searchAnalyte(100, 0.1, x, spwin = spatstat.geom::owin(c(0, 1), c(0, 1))), "every acquired")
})
