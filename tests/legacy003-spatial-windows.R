library(moleculaR)

ns <- asNamespace("moleculaR")
owin <- get("owin", ns)
as.polygonal <- get("as.polygonal", ns)
as.mask <- get("as.mask", ns)
area.owin <- get("area.owin", ns)
union.owin <- get("union.owin", ns)
setminus.owin <- get("setminus.owin", ns)
near <- function(x, y) isTRUE(all.equal(x, y, tolerance = 1e-10))

makeData <- function(coordinates, intensity, masses = 100) {
      sparseIntensityMatrix(as(Matrix::Matrix(matrix(intensity,
            nrow = nrow(coordinates)), sparse = TRUE), "dgCMatrix"),
            mzAxis = masses, coordinates = coordinates)
}
sameSignal <- function(a, b) {
      stopifnot(identical(a$n, b$n), identical(a$x, b$x), identical(a$y, b$y),
                near(a$marks$intensity, b$marks$intensity),
                identical(a$metaData$mzVals, b$metaData$mzVals),
                identical(a$metaData$mzConfirmed, b$metaData$mzConfirmed))
}

# The original rectangle reproduction, isolated from mass selection by comparison
# with the identical polygon using the same installed search implementation.
data <- makeData(data.frame(x = 1, y = 1), c(2, 10, 8), c(99.9, 100, 100.14))
rectangle <- owin(c(0.5, 1.5), c(0.5, 1.5))
polygon <- as.polygonal(rectangle)
mask <- as.mask(rectangle, dimyx = c(1, 1))
for(method in c("Gaussian", "sum", "mean", "max")) {
      a <- searchAnalyte(100, 0.1, data, spwin = rectangle, wMethod = method)
      b <- suppressWarnings(searchAnalyte(100, 0.1, data, spwin = polygon, wMethod = method))
      c <- searchAnalyte(100, 0.1, data, spwin = mask, wMethod = method)
      sameSignal(a, b)
      sameSignal(a, c)
      stopifnot(identical(a$window, rectangle), identical(b$window, polygon),
                identical(c$window, mask))
}

# A disconnected polygon with a hole retains both its domain and all detections.
main <- as.polygonal(owin(c(0.5, 15.5), c(0.5, 15.5)))
hole <- as.polygonal(owin(c(6.5, 8.5), c(6.5, 8.5)))
island <- as.polygonal(owin(c(17.5, 20.5), c(0.5, 3.5)))
complex <- union.owin(setminus.owin(main, hole), island)
coordinates <- expand.grid(x = 1:15, y = 1:15)
coordinates <- coordinates[!(coordinates$x %in% 7:8 & coordinates$y %in% 7:8), ]
coordinates <- rbind(coordinates, expand.grid(x = 18:20, y = 1:3))
data <- makeData(coordinates, seq_len(nrow(coordinates)))
result <- suppressWarnings(searchAnalyte(100, 0.1, data, spwin = complex))
stopifnot(result$n == 230, identical(result$window, complex),
          abs(area.owin(result$window) - 230) < 1e-6,
          area.owin(get("intersect.owin", ns)(result$window, hole)) < 1e-6)

# Every supported window also works for searches with no mass match or when
# confirmation excludes the analyte. Preserve the legacy empty-ppp result.
for(win in list(rectangle, polygon, mask, complex)) {
      for(confirmedOnly in c(FALSE, TRUE)) {
            empty <- suppressWarnings(searchAnalyte(if(confirmedOnly) 100 else 200,
                  0.1, data, spwin = win, verifiedMasses = 300,
                  confirmedOnly = confirmedOnly))
            stopifnot(inherits(empty, "ppp"), empty$n == 0,
                      !inherits(empty, "analytePointPattern"),
                      identical(empty$window, win),
                      identical(names(empty$marks), c("idx", "intensity")))
      }
}

# Invalid supplied windows receive an explicit class error, even with no match.
for(bad in list(NULL, "window", 1, NA_real_, list(), matrix(1), data.frame(x = 1, y = 1))) {
      for(mass in c(100, 200)) {
            error <- tryCatch(searchAnalyte(mass, 0.1, data, spwin = bad), error = identity)
            stopifnot(inherits(error, "error"), grepl("spwin.*owin", conditionMessage(error)))
      }
}

# Rectangle and polygon produce the same full scientific map with a fixed seed.
coordinates <- expand.grid(x = 1:20, y = 1:20)
intensity <- ifelse(coordinates$x %in% 7:13 & coordinates$y %in% 7:13, 10, 1)
data <- makeData(coordinates, intensity)
rectangle <- owin(c(0.5, 20.5), c(0.5, 20.5))
polygon <- as.polygonal(rectangle)
a <- searchAnalyte(100, 0.1, data, spwin = rectangle)
b <- suppressWarnings(searchAnalyte(100, 0.1, data, spwin = polygon))
sameSignal(a, b)
mapA <- suppressWarnings(probMap(a, bw = 2, seed = 1, verbose = FALSE))
mapB <- suppressWarnings(probMap(b, bw = 2, seed = 1, verbose = FALSE))
stopifnot(near(mapA$rhoMoi$v, mapB$rhoMoi$v), near(mapA$rhoCsr$v, mapB$rhoCsr$v),
          identical(mapA$hotspotIm$v > 0, mapB$hotspotIm$v > 0),
          identical(mapA$coldspotIm$v > 0, mapB$coldspotIm$v > 0))

# A representative bundled analyte agrees for equivalent rectangular domains.
example <- new.env()
data("processed-example-Data", package = "moleculaR", envir = example)
data <- createSparseMat(example$msData)
rectangle <- owin(range(data$coordinates$x) + c(-0.5, 0.5),
                   range(data$coordinates$y) + c(-0.5, 0.5))
polygon <- as.polygonal(rectangle)
width <- getFwhm(example$fwhmObj, 788.5447)
a <- searchAnalyte(788.5447, width, data, spwin = rectangle)
b <- suppressWarnings(searchAnalyte(788.5447, width, data, spwin = polygon))
sameSignal(a, b)
stopifnot(a$n == 17643)
