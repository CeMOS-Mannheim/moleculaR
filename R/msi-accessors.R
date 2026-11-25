# Internal Cardinal metadata access
#
# Cardinal remains the authoritative owner of spectra, coordinates, and runs.
# These helpers read small pixel-metadata tables; they never fetch a complete
# intensity matrix or reconstruct a global mass axis.


# Recognize both supported native representations without introducing a wrapper.
.msiIsObject <- function(x) {
    methods::is(x, "MSImagingArrays") || methods::is(x, "MSImagingExperiment")
}


# A custom project library may exist only in the parent's .libPaths(). Snow
# needs its dependencies before it can receive exported functions. Expose those
# paths during startup only, then restore the caller's environment exactly.
.msiStartBackend <- function(BPPARAM) {
    previous <- Sys.getenv("R_LIBS", unset = NA_character_)
    on.exit(if (is.na(previous)) Sys.unsetenv("R_LIBS") else Sys.setenv(R_LIBS = previous))
    Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
    BiocParallel::bpstart(BPPARAM)
}


# Return x/y in acquisition order. Row names are incidental import identifiers,
# not pixel indices. A changing z coordinate describes a volume, which must be
# reduced to one plane before using moleculaR's two-dimensional spatial methods.
.msiCoordinates <- function(x) {
    coordinates <- as.data.frame(Cardinal::coord(x))

    if (!all(c("x", "y") %in% names(coordinates))) {
        stop("MSI coordinates must contain x and y columns.")
    }

    if ("z" %in% names(coordinates) && (any(!is.finite(coordinates$z)) || length(unique(coordinates$z)) >
        1L)) {
        stop("Select a single spatial z plane before two-dimensional analysis.")
    }

    coordinates <- coordinates[, c("x", "y"), drop = FALSE]
    rownames(coordinates) <- NULL

    coordinates
}


# Use character labels so callers do not depend on factor levels left over
# after run selection. Each returned entry still corresponds to the same pixel.
.msiRun <- function(x) as.character(Cardinal::run(x))


# Metadata determine the number of acquired pixels, including analyte zeros.
.msiNumberOfPixels <- function(x) nrow(Cardinal::pixelData(x))


# Select pixels and their spectra together, preserving the requested order.
# Arrays have one index; experiments have feature rows and pixel columns.
# Cardinal 3.14.0's subsetPixels() does not handle that array distinction.
.msiSelectPixels <- function(x, pixels) {
    if (!.msiIsObject(x)) {
        stop("Expected a Cardinal MSI object.")
    }

    n <- .msiNumberOfPixels(x)

    if (is.logical(pixels)) {
        if (length(pixels) != n || anyNA(pixels)) {
            stop("Logical pixel selection must have one non-missing entry per pixel.")
        }

        pixels <- which(pixels)
    }

    if (!is.numeric(pixels) || !length(pixels) || any(!is.finite(pixels)) || any(pixels !=
        as.integer(pixels)) || any(pixels < 1 | pixels > n) || anyDuplicated(pixels)) {
        stop("Select at least one unique valid pixel index.")
    }

    if (methods::is(x, "MSImagingArrays")) {
        return(x[pixels])
    }

    x[, pixels, drop = FALSE]
}


# Validate the experiment without reading its spectral payload. Import accepts
# multiple runs and profile/unknown spectra, allowing explicit preprocessing.
# Extraction additionally requires one spatial domain, a completed processing
# queue, and known centroided status. Never guess or process implicitly here.
.msiValidate <- function(x, extraction = FALSE) {
    if (!.msiIsObject(x)) {
        stop("spData must be a Cardinal MSImagingArrays or MSImagingExperiment object.")
    }

    methods::validObject(x)

    coordinates <- .msiCoordinates(x)
    n <- .msiNumberOfPixels(x)

    if (n < 1L || nrow(coordinates) != n || !all(vapply(coordinates, is.numeric, logical(1))) ||
        any(!is.finite(as.matrix(coordinates)))) {
        stop("MSI coordinates must be finite numeric values for every acquired pixel.")
    }

    runs <- .msiRun(x)

    if (length(runs) != n || anyNA(runs) || any(!nzchar(runs))) {
        stop("Every MSI pixel must have a valid run identifier.")
    }

    # Repeated coordinates across different runs are valid; repetitions within
    # one run cannot identify a unique spatial measurement.
    if (anyDuplicated(data.frame(run = runs, coordinates))) {
        stop("Duplicate pixel coordinates within an MSI run are ambiguous.")
    }

    if (extraction) {
        if (length(unique(runs)) != 1L) {
            stop("Select a single run before extraction; use readMsi(run = ...) or native Cardinal indexing.")
        }

        if (length(ProtGenerics::processingData(x)) > 0L) {
            stop("Run Cardinal::process() to execute queued preprocessing before extraction.")
        }

        if (!isTRUE(ProtGenerics::centroided(x))) {
            stop("Extraction requires known centroided spectra; explicitly peak-pick profile data first.")
        }
    }

    invisible(x)
}
