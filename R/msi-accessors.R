# Cardinal objects stay authoritative; only pixel metadata become ordinary R data.
.msiIsObject <- function(x) {
    methods::is(x, "MSImagingArrays") || methods::is(x, "MSImagingExperiment")
}

.msiCoordinates <- function(x) {
    coordinates <- as.data.frame(Cardinal::coord(x))
    if (!all(c("x", "y") %in% names(coordinates)))
        stop("MSI coordinates must contain x and y columns.")
    if ("z" %in% names(coordinates) &&
        (any(!is.finite(coordinates$z)) || length(unique(coordinates$z)) > 1L))
        stop("Select a single spatial z plane before two-dimensional analysis.")
    coordinates <- coordinates[, c("x", "y"), drop = FALSE]
    rownames(coordinates) <- NULL
    coordinates
}

.msiRun <- function(x) as.character(Cardinal::run(x))

.msiNumberOfPixels <- function(x) nrow(Cardinal::pixelData(x))

.msiSelectPixels <- function(x, pixels) {
    if (!.msiIsObject(x)) stop("Expected a Cardinal MSI object.")
    n <- .msiNumberOfPixels(x)
    if (is.logical(pixels)) {
        if (length(pixels) != n || anyNA(pixels))
            stop("Logical pixel selection must have one non-missing entry per pixel.")
        pixels <- which(pixels)
    }
    if (!is.numeric(pixels) || !length(pixels) || any(!is.finite(pixels)) ||
        any(pixels != as.integer(pixels)) || any(pixels < 1 | pixels > n) || anyDuplicated(pixels))
        stop("Select at least one unique valid pixel index.")
    if (methods::is(x, "MSImagingArrays")) x[pixels] else x[, pixels, drop = FALSE]
}

.msiValidate <- function(x, extraction = FALSE) {
    if (!.msiIsObject(x))
        stop("spData must be a Cardinal MSImagingArrays or MSImagingExperiment object.")
    methods::validObject(x)
    coordinates <- .msiCoordinates(x)
    n <- .msiNumberOfPixels(x)
    if (n < 1L || nrow(coordinates) != n ||
        !all(vapply(coordinates, is.numeric, logical(1))) ||
        any(!is.finite(as.matrix(coordinates))))
        stop("MSI coordinates must be finite numeric values for every acquired pixel.")
    runs <- .msiRun(x)
    if (length(runs) != n || anyNA(runs) || any(!nzchar(runs)))
        stop("Every MSI pixel must have a valid run identifier.")
    if (anyDuplicated(data.frame(run = runs, coordinates)))
        stop("Duplicate pixel coordinates within an MSI run are ambiguous.")
    if (extraction) {
        if (length(unique(runs)) != 1L)
            stop("Select a single run before extraction; use readMsi(run = ...) or native Cardinal indexing.")
        if (length(ProtGenerics::processingData(x)) > 0L)
            stop("Run Cardinal::process() to execute queued preprocessing before extraction.")
        if (!isTRUE(ProtGenerics::centroided(x)))
            stop("Extraction requires known centroided spectra; explicitly peak-pick profile data first.")
    }
    invisible(x)
}
