# Run from the repository root against an installed fixed moleculaR package:
# Rscript --vanilla tests/manual/legacy001-map-comparison.R [output.csv]
# Historical extraction is loaded read-only from the verified baseline commit.
library(moleculaR)
baseline <- "a2d7f3500da5928275eb651f5c86b33beb189d67"
original <- system2("git", c("show", paste0(baseline, ":R/searchAnalyte.R")),
                    stdout = TRUE)
stopifnot(is.null(attr(original, "status")))
legacy <- new.env(parent = asNamespace("moleculaR"))
eval(parse(text = original), envir = legacy)
example <- new.env()
data("processed-example-Data", package = "moleculaR", envir = example)
spData <- createSparseMat(example$msData)
win <- createSpatialWindow(spData$coordinates, clean = FALSE, plot = FALSE)

atPixels <- function(spp) {
      index <- match(paste(spp$x, spp$y),
                     paste(spData$coordinates$x, spData$coordinates$y))
      result <- numeric(nrow(spData$coordinates))
      result[index] <- spp$marks$intensity
      result
}
classify <- function(map) {
      ifelse(map$hotspotIm$v > 0, 1L, ifelse(map$coldspotIm$v > 0, -1L, 0L))
}
results <- list()
for(mass in c(788.5447, 760.5850422936694, 496.33972664966376)) {
      width <- getFwhm(example$fwhmObj, mass)
      for(method in c("Gaussian", "sum", "mean", "max")) {
            before <- suppressWarnings(legacy$searchAnalyte(mass, width, spData,
                  wMethod = method, spwin = win))
            after <- suppressWarnings(searchAnalyte(mass, width, spData,
                  wMethod = method, spwin = win))
            oldIntensity <- atPixels(before)
            newIntensity <- atPixels(after)
            row <- data.frame(mz = mass, method = method,
                  beforePixels = before$n, afterPixels = after$n,
                  changedIntensities = sum(abs(oldIntensity - newIntensity) > 1e-12),
                  maxIntensityDelta = max(abs(oldIntensity - newIntensity)),
                  beforeTotal = sum(oldIntensity), afterTotal = sum(newIntensity))
            # A shared seed and bandwidth isolate the extraction change.
            stopifnot(before$n > 1, after$n > 1)
            oldMap <- suppressWarnings(probMap(before, bw = 2, seed = 1, verbose = FALSE))
            newMap <- suppressWarnings(probMap(after, bw = 2, seed = 1, verbose = FALSE))
            stopifnot(identical(is.na(oldMap$rhoMoi$v), is.na(newMap$rhoMoi$v)),
                      identical(before$window, after$window))
            row$maxDensityDelta <- max(abs(oldMap$rhoMoi$v - newMap$rhoMoi$v), na.rm = TRUE)
            row$maxCsrDensityDelta <- max(abs(oldMap$rhoCsr$v - newMap$rhoCsr$v), na.rm = TRUE)
            row$changedClassifications <- sum(classify(oldMap) != classify(newMap), na.rm = TRUE)
            row$beforeHot <- sum(oldMap$hotspotIm$v > 0, na.rm = TRUE)
            row$afterHot <- sum(newMap$hotspotIm$v > 0, na.rm = TRUE)
            row$beforeCold <- sum(oldMap$coldspotIm$v > 0, na.rm = TRUE)
            row$afterCold <- sum(newMap$coldspotIm$v > 0, na.rm = TRUE)
            results[[length(results) + 1L]] <- row
            cat("Compared", mass, method, "\n")
      }
}
results <- do.call(rbind, results)
print(results, row.names = FALSE, digits = 8)
args <- commandArgs(trailingOnly = TRUE)
if(length(args)) write.csv(results, args[1], row.names = FALSE)
