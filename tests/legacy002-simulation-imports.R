library(moleculaR)

stopifnot(!any(grepl("^package:spatstat", search())))
ns <- asNamespace("moleculaR")
area.owin <- get("area.owin", ns)
intersect.owin <- get("intersect.owin", ns)
disc <- get("disc", ns)
simulatespp <- get("simulatespp", ns)
# Window set operations use finite-precision polygon clipping.
near <- function(x, y) isTRUE(all.equal(x, y, tolerance = 1e-7))

# No spatstat package is attached: these calls failed in the verified baseline.
central <- createHotspotWin("centralCircle", radius = 2, backgroundLength = 20)
multiple <- createHotspotWin("multiCircle", radius = 2, backgroundLength = 20)
donut <- createHotspotWin("donut", radius = c(1, 2), backgroundLength = 20)
for(win in list(central, multiple, donut)) stopifnot(inherits(win, "owin"))

# Compare polygonal disc areas, avoiding approximation error against pi*r^2.
unitArea <- area.owin(disc(radius = 1))
stopifnot(near(area.owin(central), 4 * unitArea),
          near(area.owin(multiple), 5 * 4 * unitArea),
          near(area.owin(donut), 3 * unitArea),
          near(central$xrange, c(8, 12)), near(central$yrange, c(8, 12)),
          near(multiple$xrange, c(3, 17)), near(multiple$yrange, c(3, 17)),
          near(donut$xrange, c(8, 12)), near(donut$yrange, c(8, 12)),
          identical(central$units, multiple$units), identical(central$units, donut$units))
hole <- disc(radius = 0.5, centre = c(10, 10))
stopifnot(near(area.owin(intersect.owin(donut, hole)), 0))
reversed <- createHotspotWin("donut", radius = c(2, 1), backgroundLength = 20)
stopifnot(near(area.owin(reversed), area.owin(donut)))

# Distinct radii and overlapping circles retain their union semantics.
radii <- c(0.5, 1, 1.5, 2, 2.5)
unequal <- createHotspotWin("multiCircle", radius = radii, backgroundLength = 20)
overlap <- createHotspotWin("multiCircle", radius = 4, backgroundLength = 20)
stopifnot(near(area.owin(unequal), sum(radii^2) * unitArea),
          area.owin(overlap) < 5 * 4^2 * unitArea,
          area.owin(overlap) > 4^2 * unitArea)

# Exercise the actual simulation caller with enough intensity values to sample.
for(shape in c("centralCircle", "multiCircle", "donut")) {
      radius <- if(shape == "donut") c(1, 2) else 2
      simulate <- function() simulatespp(lambdaHotspot = 2, lambdaBackground = 0.2,
            intensityData = seq_len(10000), radius = radius, shape = shape,
            backgroundLength = 20, seed = 1, gtPlot = FALSE, sppPlot = FALSE)
      result <- simulate()
      repeatResult <- simulate()
      stopifnot(inherits(result, "analytePointPattern"), result$n > 0,
                identical(result$x, repeatResult$x), identical(result$y, repeatResult$y),
                identical(result$marks$intensity, repeatResult$marks$intensity),
                near(area.owin(result$window), 400),
                setequal(result$metaData$mzVals, c("background", "foreground")))
      foreground <- result$metaData$idx[result$metaData$mzVals == "foreground"]
      isForeground <- result$marks$idx %in% foreground
      stopifnot(any(isForeground), any(!isForeground),
                all(result$marks$intensity[isForeground] > 7500),
                all(result$marks$intensity[!isForeground] <= 7500))
}

# Repeat the simulation with measured intensities from the bundled MSI dataset.
example <- new.env()
data("processed-example-Data", package = "moleculaR", envir = example)
measuredIntensity <- unlist(lapply(example$msData[seq_len(100)], MALDIquant::intensity))
for(shape in c("centralCircle", "multiCircle", "donut")) {
      result <- simulatespp(lambdaHotspot = 2, lambdaBackground = 0.2,
            intensityData = measuredIntensity,
            radius = if(shape == "donut") c(1, 2) else 2, shape = shape,
            backgroundLength = 20, seed = 1, gtPlot = FALSE, sppPlot = FALSE)
      stopifnot(inherits(result, "analytePointPattern"), result$n > 0,
                all(result$marks$intensity %in% measuredIntensity),
                near(area.owin(result$window), 400))
}
stopifnot(!any(grepl("^package:spatstat", search())))
