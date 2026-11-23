test_that("equivalent spectra reproduce the fixed legacy extraction reference", {
    f <- legacy_fixture()
    for (x in list(fixture_cardinal(TRUE), fixture_cardinal(FALSE))) {
        for (method in names(f$patterns)) {
            actual <- searchAnalyte(f$target, f$fwhm, x, wMethod = method,
                                    verifiedMasses = c(100, 200))
            expected <- f$patterns[[method]]
            expect_s3_class(actual, "analytePointPattern")
            expect_equal(actual$x, expected$x)
            expect_equal(actual$y, expected$y)
            expect_equal(actual$marks$intensity, expected$intensity, tolerance = 1e-10)
            expect_equal(actual$metaData$mzConfirmed, expected$mzConfirmed)
            expect_equal(spatstat.geom::area.owin(actual$window), nrow(f$coordinates))
            expect_true(spatstat.geom::inside.owin(tail(f$coordinates$x, 1), tail(f$coordinates$y, 1), actual$window))
            expect_false(spatstat.geom::inside.owin(12, 3, actual$window))
        }
    }
    expect_equal(getFwhm(useDefaultFwhm(), c(249, 400, 788.5447, 1000)), f$defaultFwhm)
})

test_that("temporary legacy sparse dispatch preserves existing aggregation", {
    f <- legacy_fixture()
    sparse <- sparseIntensityMatrix(Matrix::Matrix(t(f$intensity), sparse = TRUE),
                                    f$mz, f$coordinates)
    for (method in names(f$patterns)) {
        actual <- searchAnalyte(100, 0.1, sparse, wMethod = method)
        expect_equal(actual$marks$intensity, f$patterns[[method]]$intensity, tolerance = 1e-10)
        expect_equal(actual$x, f$patterns[[method]]$x)
        expect_equal(actual$y, f$patterns[[method]]$y)
    }
})

test_that("unaligned observed-peak means differ by the documented denominator", {
    x <- fixture_cardinal(TRUE, dropZeros = TRUE)
    f <- legacy_fixture()
    expect_equal(.msiExtractIon(x, 100, 0.1, "sum")[1:5],
                 c(f$patterns$sum$intensity[1:5]))
    actual <- .msiExtractIon(x, 100, 0.1, "mean")
    count <- colSums(f$intensity[2:4, , drop = FALSE] > 0)
    expected <- numeric(ncol(f$intensity))
    expected[count > 0] <- colSums(f$intensity[2:4, , drop = FALSE])[count > 0] / count[count > 0]
    expect_equal(actual, expected)
    expect_gt(max(abs(actual[1:168] - f$patterns$mean$intensity)), 0)
})

test_that("strict mass windows quantify the nearest-index boundary difference", {
    legacy <- new.env(parent = asNamespace("moleculaR"))
    sys.source(test_path("..", "fixtures", "legacy", "searchAnalyte.R"), envir = legacy)
    mz <- c(99.9, 100, 100.14)
    x <- make_msi_arrays(list(mz), list(c(2, 10, 8)), data.frame(x = 1, y = 1))
    sparse <- sparseIntensityMatrix(Matrix::Matrix(matrix(c(2, 10, 8), 1), sparse = TRUE),
                                    mz, data.frame(x = 1, y = 1))
    win <- spatstat.geom::as.polygonal(spatstat.geom::owin(c(0.5, 1.5), c(0.5, 1.5)))
    old <- suppressWarnings(legacy$searchAnalyte(100, 0.1, sparse, spwin = win, wMethod = "sum"))
    new <- searchAnalyte(100, 0.1, x, spwin = win, wMethod = "sum")
    expect_equal(old$marks$intensity, 20)
    expect_equal(new$marks$intensity, 12)
    expect_equal(old$marks$intensity - new$marks$intensity, 8)
})

test_that("unchanged probMap reproduces seeded legacy densities and classifications", {
    f <- legacy_fixture()
    backedPath <- write_msi_pair(fixture_cardinal(TRUE), tempfile(fileext = ".imzML"))
    backed <- readMsi(backedPath, verbose = FALSE)
    expect_true(methods::is(Cardinal::intensity(backed), "matter"))
    expect_equal(Cardinal::intensity(backed)[[1]], unname(f$intensity[, 1]))
    for (x in list(fixture_cardinal(TRUE), fixture_cardinal(FALSE), backed)) {
        spp <- searchAnalyte(100, 0.1, x)
        actual <- suppressWarnings(probMap(spp, bw = 1, seed = 1, verbose = FALSE))
        expect_s3_class(actual, "molProbMap")
        expect_equal(actual$bw, f$mpm$bw)
        expect_equal(actual$rhoMoi$v, f$mpm$rhoMoi, tolerance = 1e-10)
        expect_equal(actual$rhoCsr$v, f$mpm$rhoCsr, tolerance = 1e-10)
        expect_equal(actual$hotspotIm$v, f$mpm$hotspot, tolerance = 1e-10)
        expect_equal(actual$coldspotIm$v, f$mpm$coldspot, tolerance = 1e-10)
        expect_equal(cbind(actual$hotspotpp$x, actual$hotspotpp$y), f$mpm$hotspotCoordinates)
        expect_equal(cbind(actual$coldspotpp$x, actual$coldspotpp$y), f$mpm$coldspotCoordinates)
    }
})
