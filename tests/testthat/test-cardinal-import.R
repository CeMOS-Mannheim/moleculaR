test_that("processed and continuous imzML pairs support both memory modes", {
    for (arrays in c(TRUE, FALSE)) {
        x <- if (arrays) make_msi_arrays() else make_msi_experiment()
        path <- write_msi_pair(x, tempfile(fileext = ".imzML"))
        for (memory in c(FALSE, TRUE)) {
            actual <- readMsi(path, memory = memory, verbose = FALSE)
            expect_s4_class(actual, "MSImagingArrays")
            expect_equal(.msiCoordinates(actual), .msiCoordinates(x))
            direct <- Cardinal::readMSIData(path, memory = memory, as = "MSImagingArrays", verbose = FALSE)
            expect_equal(.msiRun(actual), .msiRun(direct))
            expect_equal(Cardinal::pixelData(actual), Cardinal::pixelData(direct))
            expect_equal(Cardinal::experimentData(actual), Cardinal::experimentData(direct))
            expect_true(ProtGenerics::centroided(actual))
            expect_equal(.msiExtractIon(actual, 100, 0.1, "sum"), c(16, 3, 0))
            if (!memory) expect_true(methods::is(Cardinal::intensity(actual), "matter"))
        }
        if (!arrays) {
            shared <- readMsi(path, as = "MSImagingExperiment", verbose = FALSE)
            expect_s4_class(shared, "MSImagingExperiment")
            expect_true(methods::is(Cardinal::intensity(shared), "matter"))
            expect_equal(.msiExtractIon(shared, 100, 0.1, "sum"), c(16, 3, 0))
            expect_s4_class(readMsi(path, as = "auto", verbose = FALSE), "MSImagingExperiment")
        }
    }
})

test_that("processed profile input is imported without assuming centroid status", {
    x <- make_msi_arrays(centroided = FALSE)
    path <- write_msi_pair(x, tempfile(fileext = ".imzML"))
    actual <- readMsi(path, verbose = FALSE)
    expect_false(ProtGenerics::centroided(actual))
    expect_error(searchAnalyte(100, 0.1, actual), "centroided")
})

test_that("multi-run import supports exact run selection without mixing coordinates", {
    x <- make_msi_arrays(coordinates = data.frame(x = c(1, 1, 2), y = c(1, 1, 1)),
                         runs = c("a", "b", "a"))
    path <- write_msi_pair(x, tempfile(fileext = ".imzML"))
    all <- readMsi(path, verbose = FALSE)
    expect_equal(.msiRun(all), c("a", "b", "a"))
    expect_error(searchAnalyte(100, 0.1, all), "single run")
    selected <- readMsi(path, run = "b", verbose = FALSE)
    expect_equal(.msiCoordinates(selected), data.frame(x = 1, y = 1))
    expect_equal(.msiExtractIon(selected, 100, 0.1, "sum"), 3)
})

test_that("invalid paths, missing binaries, corrupt files and options fail", {
    expect_error(readMsi("missing.imzML"), "does not exist")
    expect_error(readMsi(NA_character_), "path must")
    expect_error(readMsi(c("a", "b")), "path must")
    bad <- tempfile(fileext = ".txt")
    writeLines("not imzML", bad)
    expect_error(readMsi(bad), "imzML files only")
    path <- tempfile(fileext = ".imzML")
    writeLines("not imzML", path)
    expect_error(readMsi(path), "paired ibd")
    writeBin(as.raw(1:10), paste0(tools::file_path_sans_ext(path), ".ibd"))
    expect_error(readMsi(path, verbose = FALSE))
    valid <- write_msi_pair(make_msi_experiment(), tempfile(fileext = ".imzML"))
    expect_error(readMsi(valid, memory = NA), "memory")
    expect_error(readMsi(valid, memory = "shared"), "memory")
    expect_error(readMsi(valid, as = "bad"), "arg")
    expect_error(readMsi(valid, run = "absent", verbose = FALSE), "not found")
    expect_error(readMsi(valid, run = NA_character_, verbose = FALSE), "run must")
    runName <- .msiRun(readMsi(valid, verbose = FALSE))[1]
    expect_equal(.msiNumberOfPixels(readMsi(valid, run = runName, verbose = FALSE)), 3)
    expect_error(readMsi(valid, parse.only = TRUE), "parse.only")
    writeBin(as.raw(1:10), paste0(tools::file_path_sans_ext(valid), ".ibd"))
    expect_error(readMsi(valid, check = TRUE, verbose = FALSE))
})
