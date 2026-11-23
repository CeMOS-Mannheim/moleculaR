.msiValidateTarget <- function(m, fwhm, method) {
    for (name in c("m", "fwhm")) {
        value <- get(name)
        if (!is.numeric(value) || length(value) != 1L ||
            !is.finite(value) || value <= 0)
            stop(name, " must be a single positive finite number.")
    }
    if (!is.character(method) || length(method) != 1L || is.na(method) ||
        !method %in% c("Gaussian", "sum", "mean", "max"))
        stop("wMethod must be one of Gaussian, sum, mean, max.")
    tolerance <- 3 * (fwhm / (2 * sqrt(2 * log(2))))
    if (!is.finite(tolerance) || tolerance <= 0)
        stop("fwhm produces an invalid search tolerance.")
    tolerance
}

.msiValidateMasses <- function(mz) {
    if (!is.numeric(mz) || any(!is.finite(mz)) || any(mz <= 0))
        stop("Measured m/z values must be positive finite numbers.")
    invisible(mz)
}

.msiAggregateIon <- function(intensity, mz, m, fwhm, method, tolerance) {
    .msiValidateMasses(mz)
    if (!is.numeric(intensity) || length(intensity) != length(mz) ||
        any(!is.finite(intensity)) || any(intensity < 0))
        stop("Retrieved intensities must be finite, non-negative, and match the m/z values.")
    selected <- which(abs(mz - m) <= tolerance)
    if (!length(selected)) return(0)
    values <- intensity[selected]
    value <- switch(method,
        Gaussian = sum(values * gaussWeight(mz[selected], m = m, fwhm = fwhm)),
        sum = sum(values), mean = mean(values), max = max(values))
    if (!is.finite(value)) stop("Aggregated analyte intensity is not finite.")
    value
}

.msiExtractIon <- function(x, m, fwhm, method = "Gaussian") {
    tolerance <- .msiValidateTarget(m, fwhm, method)
    .msiValidate(x, extraction = TRUE)
    n <- .msiNumberOfPixels(x)
    if (methods::is(x, "MSImagingExperiment")) {
        mz <- Cardinal::mz(x)
        .msiValidateMasses(mz)
        selected <- which(abs(mz - m) <= tolerance)
        if (!length(selected)) return(numeric(n))
        # Features are rows. Only the selected mass rows enter spectral iteration.
        x <- x[selected, , drop = FALSE]
    }
    values <- Cardinal::spectrapply(x, FUN = .msiAggregateIon, index = "mz",
        m = m, fwhm = fwhm, method = method, tolerance = tolerance,
        simplify = TRUE, verbose = FALSE, BPPARAM = BiocParallel::SerialParam())
    values <- as.numeric(values)
    if (length(values) != n || any(!is.finite(values)))
        stop("Cardinal extraction did not return one finite intensity per pixel.")
    values
}

.msiNearestMass <- function(m, masses, tolerance) {
    if (length(masses) == 1L && is.atomic(masses) && is.na(masses))
        return(NA_integer_)
    .msiValidateMasses(masses)
    if (!length(masses)) return(NA_integer_)
    distances <- abs(masses - m)
    nearest <- which.min(distances)
    if (distances[nearest] <= tolerance) nearest else NA_integer_
}

.searchCardinalAnalyte <- function(m, fwhm, spData, wMethod, spwin,
                                  verifiedMasses, confirmedOnly, metaData) {
    tolerance <- .msiValidateTarget(m, fwhm, wMethod)
    .msiValidate(spData, extraction = TRUE)
    if (!is.logical(confirmedOnly) || length(confirmedOnly) != 1L || is.na(confirmedOnly))
        stop("confirmedOnly must be TRUE or FALSE.")
    if (!is.list(metaData) || (length(metaData) &&
        (is.null(names(metaData)) || anyNA(names(metaData)) || any(!nzchar(names(metaData))) ||
         anyDuplicated(names(metaData)) || any(names(metaData) %in% c("idx", "mzVals", "mzConfirmed")))))
        stop("metaData must be a named list with unique names excluding idx, mzVals, mzConfirmed.")
    coordinates <- .msiCoordinates(spData)
    if (is.atomic(spwin) && length(spwin) == 1L && is.na(spwin))
        spwin <- createSpatialWindow(coordinates, plot = FALSE)
    if (!inherits(spwin, "owin")) stop("spwin must be a spatstat owin window.")
    if (!all(spatstat.geom::inside.owin(coordinates$x, coordinates$y, spwin)))
        stop("spwin must contain every acquired pixel centre.")
    spwin <- spatstat.geom::as.polygonal(spwin)
    confirmed <- !is.na(.msiNearestMass(m, verifiedMasses, tolerance))
    if (confirmedOnly && !confirmed) return(.craeteEmptySpp(spwin))
    intensities <- .msiExtractIon(spData, m, fwhm, method = wMethod)
    detected <- which(intensities > 0)
    if (!length(detected)) return(.craeteEmptySpp(spwin))
    analytePointPattern(x = coordinates$x[detected], y = coordinates$y[detected],
        win = spwin, intensity = intensities[detected], mzVals = m,
        metaData = c(list(mzConfirmed = confirmed), metaData))
}
