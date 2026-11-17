library(moleculaR)

near <- function(x, y) isTRUE(all.equal(unname(x), unname(y), tolerance = 1e-12))
fails <- function(expr) inherits(tryCatch(force(expr), error = identity), "error")

# Default FTICR calibration and preservation of the example curve's shape.
exampleData <- new.env()
data("processed-example-Data", package = "moleculaR", envir = exampleData)
reference <- exampleData$fwhmObj
model <- useDefaultFwhm()
masses <- c(150, 249, 400, 800, 1100)
stopifnot(inherits(model, "fwhm"),
          near(getFwhm(model, 249), 249 / 123000),
          near(getFwhm(model, masses) / getFwhm(reference, masses),
               rep((249 / 123000) / getFwhm(reference, 249), length(masses))),
          all(getFwhm(model, masses) > 0),
          near(model$fwhmVals, getFwhm(model, model$peaks)))

# Custom anchors, both higher and lower resolving power, and reordered names.
for(resolution in c(10000, 250000)){
      custom <- useDefaultFwhm(c(mass = 400, massRes = resolution))
      positional <- useDefaultFwhm(c(resolution, 400))
      stopifnot(near(getFwhm(custom, 400), 400 / resolution),
                near(getFwhm(custom, masses), getFwhm(positional, masses)),
                all(getFwhm(custom, masses) > 0),
                near(custom$fwhmVals, getFwhm(custom, custom$peaks)))
}

# TOF models preserve resolving power across the mass range.
for(modality in c("TOF-reflector", "TOF", "TIMS-TOF", "timsTOF")){
      resolution <- if(modality %in% c("TOF", "TOF-reflector")) 10000 else 40000
      tof <- useDefaultFwhm(modality = modality)
      custom <- useDefaultFwhm(c(30000, 600), modality = modality)
      otherAnchor <- useDefaultFwhm(c(30000, 300), modality = modality)
      fixed <- useDefaultFwhm(c(30000, 600), modality = modality, constant = TRUE)
      stopifnot(near(getFwhm(tof, masses), masses / resolution),
                near(masses / getFwhm(tof, masses), rep(resolution, length(masses))),
                near(getFwhm(custom, masses), masses / 30000),
                near(getFwhm(custom, masses), getFwhm(otherAnchor, masses)),
                near(getFwhm(fixed, masses), rep(600 / 30000, length(masses))))
}
constant <- useDefaultFwhm(c(123000, 249), constant = TRUE)
stopifnot(near(getFwhm(constant, masses), rep(249 / 123000, length(masses))))

for(bad in list(1, c(0, 249), c(-1, 249), c(123000, NA),
                c(Inf, 249), c("123000", "249"), c(R = 123000, mass = 249))){
      stopifnot(fails(useDefaultFwhm(bad)))
}
stopifnot(fails(useDefaultFwhm(modality = "unknown")),
          fails(useDefaultFwhm(modality = NA_character_)),
          fails(useDefaultFwhm(constant = NA)),
          fails(useDefaultFwhm(constant = 1)))

# The model can be serialized without relying on a loaded example dataset.
path <- tempfile(fileext = ".rds")
saveRDS(model, path)
restored <- readRDS(path)
unlink(path)
stopifnot(near(getFwhm(restored, masses), getFwhm(model, masses)),
          !exists("fwhmObj", envir = .GlobalEnv, inherits = FALSE),
          !exists("msData", envir = .GlobalEnv, inherits = FALSE))
