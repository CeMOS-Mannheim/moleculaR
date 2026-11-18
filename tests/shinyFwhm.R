if(requireNamespace("testthat", quietly = TRUE)) local({
      library(moleculaR)
      library(shiny)
      appDir <- system.file("shiny", package = "moleculaR")
      oldDir <- setwd(appDir)
      on.exit(setwd(oldDir))
      app <- new.env(parent = globalenv())
      expressions <- parse(file.path(appDir, "package-app.R"))
      serverDefinition <- Filter(function(e) {
            is.call(e) && identical(e[[1]], as.name("<-")) &&
                  identical(e[[2]], as.name("server"))
      }, expressions)
      eval(serverDefinition[[1]], envir = app)

      # Exercise upload preparation and model selection without a large imzML fixture.
      coordinates <- expand.grid(x = 1:3, y = 1:3)
      spectra <- lapply(seq_len(nrow(coordinates)), function(i) {
            MALDIquant::createMassPeaks(mass = c(300, 400, 500), intensity = c(1, 2, 3),
                  metaData = list(imaging = list(pos = unlist(coordinates[i, ]))))
      })
      reads <- 0L
      models <- list()
      estimate <- moleculaR::estimateFwhm
      testthat::local_mocked_bindings(readCentrData = function(path, ...) {
            stopifnot(file.exists(path), file.exists(sub("imzML$", "ibd", path)))
            reads <<- reads + 1L
            spectra
      }, estimateFwhm = function(...) {
            model <- estimate(...)
            models[[length(models) + 1L]] <<- model
            model
      }, .package = "moleculaR")
      uploadPaths <- c(tempfile(), tempfile())
      file.create(uploadPaths)
      on.exit(unlink(uploadPaths), add = TRUE)
      upload <- data.frame(name = c("sample.IMZML", "sample.IBD"), datapath = uploadPaths,
                            size = 0, type = "application/octet-stream")
      profile <- system.file("extdata", "pos-XIII-82492-singleSpectrum.tsv", package = "moleculaR")

      shiny::testServer(app$server, {
            session$setInputs(fwhmSource = "default", fwhmModality = "TIMS-TOF",
                              customMassRes = FALSE, imzmlFile = upload, loadData = 1)
            stopifnot(reads == 1L,
                      isTRUE(all.equal(moleculaR::getFwhm(models[[1]], 800), 800 / 40000)),
                      grepl("TIMS-TOF", output$message), grepl("40000", output$message))

            # Custom calibration and repeated loading keep the original uploads intact.
            session$setInputs(customMassRes = TRUE, fwhmMassRes = 30000,
                              fwhmMass = 600, loadData = 2)
            stopifnot(reads == 2L, all(file.exists(uploadPaths)),
                      isTRUE(all.equal(moleculaR::getFwhm(models[[2]], 600), 600 / 30000)))

            # The measured route requires a spectrum and displays a useful message.
            session$setInputs(fwhmSource = "spectrum", loadData = 3)
            stopifnot(reads == 2L, grepl("upload a spectrum", output$message))
            session$setInputs(spectrFile = data.frame(name = "profile.tsv", datapath = profile),
                              loadData = 4)
            stopifnot(reads == 3L, grepl("estimated from the uploaded spectrum", output$message))

            # Switching back explicitly uses defaults even if a spectrum remains uploaded.
            session$setInputs(fwhmSource = "default", customMassRes = FALSE,
                              fwhmModality = "TOF-reflector", loadData = 5)
            stopifnot(reads == 4L,
                      isTRUE(all.equal(moleculaR::getFwhm(models[[4]], 800), 800 / 10000)),
                      grepl("approximate", output$message))
      })
      # A fresh session cannot use another session's initialized model.
      shiny::testServer(app$server, {
            session$setInputs(renderMPM = 1)
            stopifnot(grepl("No FWHM data", output$message))
      })
})
