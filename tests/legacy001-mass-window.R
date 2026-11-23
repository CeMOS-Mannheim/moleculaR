library(moleculaR)

# Resolve geometry through moleculaR's dynamic spatstat imports.
owin <- get("owin", asNamespace("moleculaR"))
as.polygonal <- get("as.polygonal", asNamespace("moleculaR"))
area.owin <- get("area.owin", asNamespace("moleculaR"))
near <- function(x, y) isTRUE(all.equal(x, y, tolerance = 1e-12))
methods <- c("Gaussian", "sum", "mean", "max")

fixture <- function(masses, intensity) {
      intensity <- as.matrix(intensity)
      sparseIntensityMatrix(
            spmat = as(Matrix::Matrix(intensity, sparse = TRUE), "dgCMatrix"),
            mzAxis = masses,
            coordinates = data.frame(x = seq_len(nrow(intensity)), y = 1))
}
windowFor <- function(data) {
      as.polygonal(owin(range(data$coordinates$x) + c(-0.5, 0.5),
                        range(data$coordinates$y) + c(-0.5, 0.5)))
}
checkExtraction <- function(data, mass, width, expected) {
      win <- windowFor(data)
      for(method in methods) {
            result <- suppressWarnings(searchAnalyte(mass, width, data,
                  wMethod = method, spwin = win))
            positive <- which(expected[[method]] > 0)
            stopifnot(inherits(result, "ppp"), result$n == length(positive),
                      near(result$x, data$coordinates$x[positive]),
                      near(result$marks$intensity, expected[[method]][positive]),
                      identical(result$window, win))
      }
}

# Historical reproduction: legacy sum was 20, including the outside 100.14 peak.
data <- fixture(c(99.90, 100, 100.14), matrix(c(2, 10, 8), nrow = 1))
sigma <- 0.1 / (2 * sqrt(2 * log(2)))
checkExtraction(data, 100, 0.1, list(
      Gaussian = 10 + 2 * exp(-0.5 * (0.1 / sigma)^2),
      sum = 12, mean = 6, max = 10))

# Exactly representable three-sigma endpoints and their immediate neighbours.
# Outside peaks are deliberately intense, including a pixel with ONLY outside
# peaks. The fourth acquired pixel is all-zero and must remain in the domain.
width <- 2 * sqrt(2 * log(2))
epsilon <- 32 * .Machine$double.eps
data <- fixture(c(13 - epsilon, 13, 15.9, 16, 18.1, 19, 19 + epsilon),
      rbind(c(1e12, 2, 3, 5, 7, 11, 1e12),
            c(0, 0, 0, 13, 0, 0, 0),
            c(1e12, 0, 0, 0, 0, 0, 1e12),
            rep(0, 7)))
checkExtraction(data, 16, width, list(
      Gaussian = c(13 * exp(-4.5) + 3 * exp(-0.5 * 0.1^2) +
                         5 + 7 * exp(-0.5 * 2.1^2), 13, 0, 0),
      sum = c(28, 13, 0, 0), mean = c(28 / 5, 13 / 5, 0, 0),
      max = c(11, 13, 0, 0)))

# A coarse/nonuniform axis must not pull in the nearest outside neighbours.
data <- fixture(c(1, 12.9, 15.75, 16.5, 19.1, 100),
                matrix(c(1e12, 1e12, 4, 8, 1e12, 1e12), nrow = 1))
checkExtraction(data, 16, width, list(
      Gaussian = 4 * exp(-0.5 * 0.25^2) + 8 * exp(-0.5 * 0.5^2),
      sum = 12, mean = 6, max = 8))

# One selected feature must retain matrix dimensions for every aggregation.
data <- fixture(c(1, 16, 100), rbind(c(1e12, 9, 1e12), c(0, 0, 0)))
checkExtraction(data, 16, width, setNames(rep(list(c(9, 0)), 4), methods))

# No columns selected, an empty feature axis, or zero intensities all retain the
# legacy empty-ppp convention and its marks and window.
for(data in list(
      fixture(c(12.9, 19.1), matrix(c(1e12, 1e12), nrow = 1)),
      fixture(numeric(), matrix(numeric(), nrow = 1, ncol = 0)),
      fixture(c(15, 16, 17), matrix(0, nrow = 1, ncol = 3)))) {
      checkExtraction(data, 16, width, setNames(rep(list(0), 4), methods))
      empty <- suppressWarnings(searchAnalyte(16, width, data, spwin = windowFor(data)))
      stopifnot(!inherits(empty, "analytePointPattern"),
                identical(names(empty$marks), c("idx", "intensity")))
}

# Auto-generated tissue windows include acquired pixels with no selected signal.
data <- fixture(c(16, 20), rbind(c(5, 0), c(0, 1e12), c(0, 0), c(3, 0)))
data$coordinates <- expand.grid(x = 1:2, y = 1:2)
result <- searchAnalyte(16, width, data, wMethod = "sum")
stopifnot(result$n == 2, abs(area.owin(result$window) - 4) < 1e-7)

# Confirmation, metadata and batch extraction retain the corrected intensities.
win <- windowFor(data)
confirmed <- suppressWarnings(searchAnalyte(16, width, data, spwin = win,
      verifiedMasses = 16, confirmedOnly = TRUE, metaData = list(label = "fixture")))
stopifnot(confirmed$metaData$mzConfirmed, confirmed$metaData$label == "fixture")
unconfirmed <- suppressWarnings(searchAnalyte(16, width, data, spwin = win,
      verifiedMasses = 100, confirmedOnly = TRUE))
stopifnot(unconfirmed$n == 0)
batch <- suppressWarnings(batchAnalyteSearch(data,
      fwhm(function(m) rep(width, length(m))), spwin = win,
      m = 16, wMethod = "sum", verbose = FALSE))
stopifnot(near(batch$marks$intensity, c(5, 3)))
