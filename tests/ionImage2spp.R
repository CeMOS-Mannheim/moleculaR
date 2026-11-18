library(moleculaR)

fails <- function(expr) inherits(tryCatch(force(expr), error = identity), "error")
same <- function(x, y) isTRUE(all.equal(x, y, check.attributes = FALSE))

# Asymmetric values and dimensions detect accidental transposes or vertical flips.
image <- matrix(c(1, 0, NA, 2, 3, 4, NA, 5, 6, 7, NA, 8), nrow = 3)
spp <- ionImage2spp(image, mzVal = 788.5)
idx <- which(!is.na(image) & image > 0, arr.ind = TRUE)
stopifnot(inherits(spp, "analytePointPattern"), inherits(spp, "ppp"),
          identical(spp$window$type, "polygonal"),
          same(spp$x, idx[, 1]), same(spp$y, idx[, 2]),
          same(spp$marks$intensity, image[idx]), same(spp$metaData$mzVals, 788.5),
          same(spatstat.geom::area.owin(spp$window), sum(!is.na(image))))
withZeros <- ionImage2spp(image, zero.rm = FALSE)
stopifnot(withZeros$n == sum(!is.na(image)), any(withZeros$marks$intensity == 0),
          same(spatstat.geom::area.owin(withZeros$window), spatstat.geom::area.owin(spp$window)))
transposed <- ionImage2spp(t(image), orientation = "yx", mzVal = 788.5)
key <- function(p) paste(p$x, p$y, p$marks$intensity, sep = ":")
stopifnot(identical(sort(key(spp)), sort(key(transposed))))

# Missing pixels define holes and disconnected tissue components.
mask <- matrix(NA_real_, 8, 5)
mask[1:3, 1:3] <- 1
mask[2, 2] <- NA
mask[7:8, 4:5] <- 2
pattern <- ionImage2spp(mask)
stopifnot(same(spatstat.geom::area.owin(pattern$window), 12),
          !spatstat.geom::inside.owin(2, 2, pattern$window),
          !spatstat.geom::inside.owin(5, 3, pattern$window),
          spatstat.geom::inside.owin(8, 5, pattern$window))
commonWin <- createSpatialWindow(expand.grid(x = 1:8, y = 1:5), plot = FALSE)
shared <- ionImage2spp(mask, win = commonWin)
stopifnot(same(spatstat.geom::area.owin(shared$window), 40))
rectangular <- createSpatialWindow(data.frame(x = c(10, 12), y = c(20, 23)),
                                    pixelSize = c(2, 3), plot = FALSE)
stopifnot(same(spatstat.geom::area.owin(rectangular), 12),
          fails(createSpatialWindow(data.frame(x = 1:2, y = 1:2), pixelSize = 0, plot = FALSE)),
          fails(createSpatialWindow(data.frame(x = c(1, 1.5), y = 1:2), pixelSize = 1, plot = FALSE)))
stopifnot(fails(ionImage2spp(mask, win = "window")),
          fails(ionImage2spp(mask, win = spatstat.geom::owin(c(0, 3), c(0, 3)))))

# Singleton dimensions are retained; multiple mass layers require selection.
one <- array(image, c(3, 4, 1))
attr(one, "center") <- 788.5
stopifnot(identical(key(ionImage2spp(one)), key(spp)),
          same(ionImage2spp(one)$metaData$mzVals, 788.5))
layers <- array(c(image, image * 2), c(3, 4, 2))
attr(layers, "center") <- c(788.5, 810.6)
second <- ionImage2spp(layers, slice = 2)
stopifnot(fails(ionImage2spp(layers)),
          same(second$marks$intensity, spp$marks$intensity * 2),
          same(second$metaData$mzVals, 810.6),
          same(ionImage2spp(layers, slice = 2, mzVal = 900)$metaData$mzVals, 900))
for(bad in list(0, 3, 1.5, NA, Inf, "1", c(1, 2))) {
      stopifnot(fails(ionImage2spp(layers, slice = bad)))
}
stopifnot(fails(ionImage2spp(image, slice = 1)),
          ionImage2spp(matrix(1:4, 1))$n == 4,
          ionImage2spp(array(1:4, c(4, 1, 1)))$n == 4)
isolated <- matrix(NA_real_, 5, 5)
isolated[cbind(c(1, 3, 5), c(1, 3, 5))] <- 1:3
stopifnot(same(spatstat.geom::area.owin(ionImage2spp(isolated)$window), 3),
          same(spatstat.geom::area.owin(ionImage2spp(matrix(1))$window), 1))

# Real msiSlices output: preserve [x, y, mass] coordinates and mass metadata.
coords <- expand.grid(x = 2:4, y = 4:8)
spectra <- lapply(seq_len(nrow(coords)), function(i) {
      MALDIquant::createMassPeaks(mass = c(400, 800), intensity = c(i, i * 2),
            metaData = list(imaging = list(pos = unlist(coords[i, ]))))
})
ionArray <- MALDIquant::msiSlices(spectra, center = c(400, 800), tolerance = 0.01, adjust = FALSE)
actual <- ionImage2spp(ionArray, slice = 2)
stopifnot(same(actual$x, coords$x), same(actual$y, coords$y),
          same(actual$marks$intensity, seq_len(nrow(coords)) * 2),
          same(actual$metaData$mzVals, 800))

# Bad intensities and shapes fail before spatial analysis.
for(bad in list(matrix(NA_real_, 2, 2), matrix(0, 2, 2), matrix(-1, 2, 2),
                matrix(Inf, 2, 2), matrix("1", 2, 2), matrix(TRUE, 2, 2),
                1:4, array(1, c(2, 2, 1, 1)), matrix(numeric(), 0, 2))) {
      stopifnot(fails(ionImage2spp(bad)))
}
stopifnot(fails(ionImage2spp(image, orientation = "invalid")),
          fails(ionImage2spp(image, zero.rm = NA)),
          fails(ionImage2spp(image, mzVal = c(400, 800))))

# End-to-end mapping uses the original image intensities and tissue domain.
image <- matrix(NA_real_, 20, 15)
image[3:18, 3:13] <- 1
image[7:13, 6:10] <- 10
spp <- ionImage2spp(image, mzVal = 788.5)
mpm <- probMap(spp, bw = 1, seed = 1, verbose = FALSE)
stopifnot(inherits(mpm, "molProbMap"),
          same(mpm$sppMoi$marks$intensity, spp$marks$intensity),
          all(is.finite(mpm$rhoMoi$v[!is.na(mpm$rhoMoi$v)])))
