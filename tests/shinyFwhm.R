# Installed-app smoke tests complement the detailed source/native server suite.
if (requireNamespace("testthat", quietly = TRUE)) local({
    library(moleculaR)
    library(shiny)
    appDir <- system.file("shiny", package = "moleculaR")
    example <- system.file("extdata", "shiny-example", "example.imzML", package = "moleculaR")
    original <- tools::md5sum(c(example, sub("imzML$", "ibd", example)))
    for (name in c("package-app", "example-app")) {
        app <- new.env(parent = globalenv())
        sys.source(file.path(appDir, paste0(name, ".R")), envir = app)
        root <- NULL
        shiny::testServer(app$server, {
            session$flushReact()
            root <<- modules$example$state$root
            session$setInputs(`example-fwhmSource` = "default", `example-fwhmModality` = "TIMS-TOF",
                `example-normalization` = "none", `example-loadData` = 1)
            state <- modules$example$state
            stopifnot(state$ready, methods::is(state$msData, "MSImagingArrays"),
                methods::is(Cardinal::intensity(state$msData), "matter"),
                isTRUE(all.equal(getFwhm(state$fwhmObj, 800), 800 / 40000)))
            session$setInputs(`example-mz` = 788.5447, `example-wMethod` = "Gaussian",
                `example-confirmedOnly` = TRUE, `example-bandwidthMode` = "fixed", `example-bandwidth` = 2,
                `example-seed` = 1)
            session$setInputs(`example-renderMPM` = 1)
            stopifnot(!is.null(state$result), state$result$pointPattern$n > 100,
                state$result$pointPattern$metaData$mzVals == 788.5447)
            if (name == "package-app") stopifnot(!modules$main$state$ready, is.null(modules$main$state$result),
                !identical(modules$main$state$root, root))
        })
        stopifnot(!dir.exists(root), identical(tools::md5sum(c(example, sub("imzML$", "ibd", example))), original))
    }
})
