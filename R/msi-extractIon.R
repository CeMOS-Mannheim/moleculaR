# Internal single-analyte extraction
#
# Data flow: validate the query -> extract one intensity per acquired pixel ->
# select positive detections -> construct the existing spatial point pattern.
# Scientific choices, including mass boundaries and mean denominators, are
# recorded in docs/Cardinal-migration/scientific-decisions.md (D002-D005).


# Validate one query and return its absolute mass tolerance, in Da. Keeping this
# conversion in one place ensures extraction and annotation matching agree.
.msiValidateTarget <- function(m, fwhm, method) {
    for (name in c("m", "fwhm")) {
        value <- get(name)
        if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value <= 0) {
            stop(name, " must be a single positive finite number.")
        }
    }

    if (!is.character(method) || length(method) != 1L || is.na(method) || !method %in% c("Gaussian",
        "sum", "mean", "max")) {
        stop("wMethod must be one of Gaussian, sum, mean, max.")
    }

    # FWHM = 2 * sqrt(2 * log(2)) * sigma; search three sigma on each side.
    tolerance <- 3 * (fwhm / (2 * sqrt(2 * log(2))))

    if (!is.finite(tolerance) || tolerance <= 0) {
        stop("fwhm produces an invalid search tolerance.")
    }

    tolerance
}


# Empty spectra are valid. Duplicate and unsorted measured masses are also
# valid: the extraction rule evaluates each observation, without alignment.
.msiValidateMasses <- function(mz) {
    if (!is.numeric(mz) || any(!is.finite(mz)) || any(mz <= 0)) {
        stop("Measured m/z values must be positive finite numbers.")
    }

    invisible(mz)
}


# Cardinal calls this worker once per spectrum. Return one finite number,
# including zero for an empty interval, so the result retains pixel order.
.msiAggregateIon <- function(intensity, mz, m, fwhm, method, tolerance) {
    .msiValidateMasses(mz)

    if (!is.numeric(intensity) || length(intensity) != length(mz) || any(!is.finite(intensity)) ||
        any(intensity < 0)) {
        stop("Retrieved intensities must be finite, non-negative, and match the m/z values.")
    }

    # Endpoints are inclusive. Unlike the legacy nearest-index interval, this
    # never extends the window to include a peak outside the numerical bound.
    selected <- which(abs(mz - m) <= tolerance)

    if (!length(selected)) {
        return(0)
    }

    values <- intensity[selected]

    # For arrays, the denominator of mean is the observed peak count. Shared
    # experiments supply every selected feature, including zero-filled gaps.
    value <- switch(
        method,
        Gaussian = sum(values * gaussWeight(mz[selected], m = m, fwhm = fwhm)),
        sum = sum(values),
        mean = mean(values),
        max = max(values)
    )

    if (!is.finite(value)) {
        stop("Aggregated analyte intensity is not finite.")
    }

    value
}


# Extract a numeric vector with exactly one entry per acquired pixel. Cardinal
# manages file access and chunks; only selected feature rows are needed for a
# shared-axis experiment. Individual-axis arrays are visited spectrum by spectrum.
.msiExtractIon <- function(x, m, fwhm, method = "Gaussian") {
    tolerance <- .msiValidateTarget(m, fwhm, method)
    .msiValidate(x, extraction = TRUE)
    n <- .msiNumberOfPixels(x)

    if (methods::is(x, "MSImagingExperiment")) {
        mz <- Cardinal::mz(x)
        .msiValidateMasses(mz)
        selected <- which(abs(mz - m) <= tolerance)

        if (!length(selected)) {
            return(numeric(n))
        }

        # Features are rows. Only the selected mass rows enter spectral iteration.
        x <- x[selected, , drop = FALSE]
    }

    # Serial execution is intentional for this milestone. Leave chunk sizes to
    # Cardinal's configuration; parallel/chunk tuning needs separate benchmarks.
    values <- Cardinal::spectrapply(
        x,
        FUN = .msiAggregateIon,
        index = "mz",
        m = m,
        fwhm = fwhm,
        method = method,
        tolerance = tolerance,
        simplify = TRUE,
        verbose = FALSE,
        BPPARAM = BiocParallel::SerialParam()
    )

    values <- as.numeric(values)

    if (length(values) != n || any(!is.finite(values))) {
        stop("Cardinal extraction did not return one finite intensity per pixel.")
    }

    values
}


# Return the original index of the nearest verified mass within tolerance.
# which.min() breaks ties by input order. The public NA sentinel and an empty
# annotation vector both mean that there is no confirmation for the target.
.msiNearestMass <- function(m, masses, tolerance) {
    if (length(masses) == 1L && is.atomic(masses) && is.na(masses)) {
        return(NA_integer_)
    }

    .msiValidateMasses(masses)

    if (!length(masses)) {
        return(NA_integer_)
    }

    distances <- abs(masses - m)
    nearest <- which.min(distances)

    if (distances[nearest] <= tolerance) {
        return(nearest)
    }

    NA_integer_
}


# Adapt Cardinal extraction to moleculaR's existing public spatial interface.
# No Cardinal object enters probMap(): its input remains analytePointPattern.
.searchCardinalAnalyte <- function(m, fwhm, spData, wMethod, spwin,
                                   verifiedMasses, confirmedOnly, metaData) {
    tolerance <- .msiValidateTarget(m, fwhm, wMethod)
    .msiValidate(spData, extraction = TRUE)

    if (!is.logical(confirmedOnly) || length(confirmedOnly) != 1L || is.na(confirmedOnly)) {
        stop("confirmedOnly must be TRUE or FALSE.")
    }

    # Reserved names link point marks to analyte metadata. Caller-supplied names
    # must not shadow those identifiers or the annotation-confirmation field.
    invalidMetadata <- !is.list(metaData)

    if (!invalidMetadata && length(metaData)) {
        metadataNames <- names(metaData)
        reservedNames <- c("idx", "mzVals", "mzConfirmed")

        invalidMetadata <- is.null(metadataNames) ||
            anyNA(metadataNames) ||
            any(!nzchar(metadataNames)) ||
            anyDuplicated(metadataNames) ||
            any(metadataNames %in% reservedNames)
    }

    if (invalidMetadata) {
        stop("metaData must be a named list with unique names excluding idx, mzVals, mzConfirmed.")
    }

    # Build/validate the window before discarding any analyte zeros. Acquired
    # zero pixels belong to the tissue domain even though they are not detections.
    coordinates <- .msiCoordinates(spData)

    if (is.atomic(spwin) && length(spwin) == 1L && is.na(spwin)) {
        spwin <- createSpatialWindow(coordinates, plot = FALSE)
    }

    if (!inherits(spwin, "owin")) {
        stop("spwin must be a spatstat owin window.")
    }

    if (!all(spatstat.geom::inside.owin(coordinates$x, coordinates$y, spwin))) {
        stop("spwin must contain every acquired pixel centre.")
    }

    spwin <- spatstat.geom::as.polygonal(spwin)
    confirmed <- !is.na(.msiNearestMass(m, verifiedMasses, tolerance))

    # An unconfirmed query can return before reading any spectral values.
    if (confirmedOnly && !confirmed) {
        return(.craeteEmptySpp(spwin))
    }

    intensities <- .msiExtractIon(spData, m, fwhm, method = wMethod)
    detected <- which(intensities > 0)

    if (!length(detected)) {
        return(.craeteEmptySpp(spwin))
    }

    analytePointPattern(
        x = coordinates$x[detected],
        y = coordinates$y[detected],
        win = spwin,
        intensity = intensities[detected],
        mzVals = m,
        metaData = c(list(mzConfirmed = confirmed), metaData)
    )
}
