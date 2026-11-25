# Shared native preprocessing utilities. No global spectral conversion occurs.

.msiCheckProcessing <- function(x) {
    .msiValidate(x)
    if (length(ProtGenerics::processingData(x))) {
        stop("Run Cardinal::process() to execute queued preprocessing first.")
    }
    invisible(x)
}

.numericSpectrum <- function(mz, intensity) {
    if (!is.numeric(mz) || !is.numeric(intensity) || !is.null(dim(mz)) || !is.null(dim(intensity)) || !length(mz) ||
        length(mz) != length(intensity) || any(!is.finite(mz)) || any(mz <= 0) ||
        any(!is.finite(intensity)) || any(intensity < 0)) {
        stop("A profile requires equal non-empty numeric mz/intensity vectors with positive finite masses and finite non-negative intensities.")
    }
    if (is.unsorted(mz, strictly = TRUE)) {
        stop("Profile m/z values must be strictly increasing, without duplicates.")
    }
    structure(list(mz = as.numeric(mz), intensity = as.numeric(intensity)),
              class = "moleculaRSpectrum")
}

.msiSpectrumVectors <- function(x, pixel) {
    if (methods::is(x, "MSImagingArrays")) {
        list(mz = as.numeric(Cardinal::mz(x)[[pixel]]),
             intensity = as.numeric(Cardinal::intensity(x)[[pixel]]))
    } else {
        list(mz = as.numeric(Cardinal::mz(x)),
             intensity = as.numeric(Cardinal::intensity(x)[, pixel]))
    }
}

# Reference medians and frequencies need feature rows across pixels. Bound their
# dense workspace by Cardinal's byte setting, or a 16 MiB target, not all rows.
.msiFeatureBlocks <- function(x, featureChunkSize = NULL) {
    n <- length(Cardinal::mz(x))
    if (is.null(featureChunkSize)) {
        budget <- Cardinal::getCardinalChunksize()
        if (!inherits(budget, "size_bytes") || length(budget) != 1L ||
            !is.finite(budget) || budget <= 0) budget <- 16 * 1024^2
        budget <- as.numeric(budget)
        featureChunkSize <- max(1, floor(budget / (8 * .msiNumberOfPixels(x))))
    }
    if (!is.numeric(featureChunkSize) || length(featureChunkSize) != 1L ||
        !is.finite(featureChunkSize) || featureChunkSize < 1 ||
        featureChunkSize != floor(featureChunkSize)) {
        stop("featureChunkSize must be a positive integer number of feature rows.")
    }
    if (!n) return(list())
    starts <- seq.int(1, n, by = featureChunkSize)
    lapply(starts, function(start) seq.int(start, min(n, start + featureChunkSize - 1)))
}

.msiValidateIntensities <- function(intensity, mz = NULL) {
    if (!is.numeric(intensity) || any(!is.finite(intensity)) || any(intensity < 0) ||
        (!is.null(mz) && length(intensity) != length(mz))) {
        stop("Retrieved intensities must be finite, non-negative, and match the m/z values.")
    }
    invisible(intensity)
}

.msiRecordProcessing <- function(x, label, details) {
    processing <- ProtGenerics::processingData(x)
    attr(processing, "processed") <- c(attr(processing, "processed"),
                                      setNames(list(details), label))
    ProtGenerics::processingData(x) <- processing
    methods::validObject(x)
    x
}

.msiFoldReference <- function(x, meanMethod, zeroAsMissing, featureChunkSize) {
    if (!methods::is(x, "MSImagingExperiment")) {
        stop("medianFoldChange requires an aligned MSImagingExperiment; align arrays explicitly first.")
    }
    if (!isTRUE(ProtGenerics::centroided(x))) {
        stop("medianFoldChange requires known centroided shared features.")
    }
    meanMethod <- match.arg(meanMethod, c("mean", "median", "sum"))
    reference <- numeric(length(Cardinal::mz(x)))
    for (rows in .msiFeatureBlocks(x, featureChunkSize)) {
        block <- as.matrix(Cardinal::intensity(x)[rows, , drop = FALSE])
        .msiValidateIntensities(block)
        reference[rows] <- apply(block, 1, function(values) {
            if (zeroAsMissing) values <- values[values > 0]
            if (!length(values)) return(0)
            switch(meanMethod, mean = mean(values), median = stats::median(values), sum = sum(values))
        })
    }
    if (any(!is.finite(reference))) stop("Fold-change reference contains a non-finite value.")
    reference
}

# Scalar factor plus whether the acquired spectrum has positive signal. Raw
# missing/zero standard factors are retained for provenance, not imputed.
.msiNormalizationFactor <- function(intensity, mz, method, reference, standard, tolerance,
                                     zeroAsMissing, scale) {
    .msiValidateMasses(mz)
    .msiValidateIntensities(intensity, mz)
    positive <- any(intensity > 0)
    if (!positive && method != "IS") return(c(factor = 1, positive = 0))
    factor <- switch(method,
        TIC = sum(intensity),
        RMS = sqrt(mean(intensity^2)),
        median = stats::median(intensity),
        IS = {
            if (is.unsorted(mz)) stop("Reference normalization requires sorted m/z values.")
            index <- if (length(mz)) matter::bsearch(standard, mz, nearest = TRUE) else NA_integer_
            if (is.na(index) || abs(mz[index] - standard) > tolerance) NA_real_ else intensity[index]
        },
        medianFoldChange = {
            keep <- reference > 0
            if (zeroAsMissing) keep <- keep & intensity > 0
            if (!any(keep)) NA_real_ else stats::median(intensity[keep] / reference[keep])
        })
    if (method != "IS" && positive && (!is.finite(factor) || factor <= 0)) {
        stop("Normalization factor is undefined or non-positive for a nonzero spectrum.")
    }
    c(factor = if (positive) unname(factor / scale) else 1, positive = as.numeric(positive))
}

.msiNormalizeSpectrum <- function(intensity, mz, method, reference, standard, tolerance,
                                 zeroAsMissing, scale) {
    summary <- .msiNormalizationFactor(intensity, mz, method, reference, standard,
                                        tolerance, zeroAsMissing, scale)
    if (!summary["positive"]) return(intensity)
    value <- switch(method,
        TIC = matter::rescale_sum(intensity, scale = scale),
        RMS = matter::rescale_rms(intensity, scale = scale),
        IS = {
            if (is.na(summary["factor"]) || summary["factor"] == 0) {
                rep.int(0, length(intensity))
            } else {
                matter::rescale_ref(intensity, ref = standard, scale = scale, domain = mz)
            }
        },
        intensity / summary["factor"])
    .msiValidateIntensities(value, mz)
    as.numeric(value)
}

.normalizeCardinal <- function(x, method, meanMethod, mz, tol, scale, zeroAsMissing,
                               outfile, featureChunkSize, BPPARAM, verbose) {
    .msiCheckProcessing(x)
    method <- match.arg(method, c("TIC", "RMS", "median", "IS", "medianFoldChange"))
    if (!is.logical(zeroAsMissing) || length(zeroAsMissing) != 1L || is.na(zeroAsMissing)) {
        stop("zeroAsMissing must be TRUE or FALSE.")
    }
    if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) stop("verbose must be TRUE or FALSE.")
    if (is.logical(scale) && identical(scale, NA)) scale <- NA_real_
    if (!is.numeric(scale) || length(scale) != 1L || is.nan(scale) ||
        (!is.na(scale) && (!is.finite(scale) || scale <= 0))) stop("scale must be NA or a positive finite number.")
    if (is.na(scale)) {
        scale <- if (method %in% c("TIC", "RMS")) {
            if (methods::is(x, "MSImagingArrays")) max(lengths(Cardinal::mz(x))) else length(Cardinal::mz(x))
        } else 1
        scale <- max(1, scale)
    }
    if (method == "IS") {
        if (!is.numeric(mz) || length(mz) != 1L || !is.finite(mz) || mz <= 0 ||
            !is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol < 0) {
            stop("IS requires a positive finite mz and a non-negative finite absolute tol.")
        }
    }
    if (method == "medianFoldChange" && length(unique(.msiRun(x))) != 1L) {
        stop("Select one run before estimating a fold-change reference.")
    }
    reference <- if (method == "medianFoldChange")
        .msiFoldReference(x, meanMethod, zeroAsMissing, featureChunkSize) else numeric()
    backend <- .msiBatchBackend(1, BPPARAM)
    started <- !BiocParallel::bpisup(backend)
    if (started) {
        .msiStartBackend(backend)
        on.exit(BiocParallel::bpstop(backend), add = TRUE)
    }
    factors <- Cardinal::spectrapply(x, .msiNormalizationFactor, index = "mz",
        method = method, reference = reference, standard = mz, tolerance = tol,
        zeroAsMissing = zeroAsMissing, scale = scale, simplify = TRUE,
        BPPARAM = backend, verbose = FALSE)
    factors <- matrix(as.numeric(factors), nrow = 2)

    # Direct native backing avoids Cardinal 3.14/matter recursion when process()
    # combines thousands of separately backed results. Only one retrieved
    # spectrum is resident during this write pass; original m/z storage is kept.
    if (is.null(outfile)) outfile <- tempfile("moleculaR-normalized-", fileext = ".bin")
    if (!is.character(outfile) || length(outfile) != 1L || is.na(outfile) || !nzchar(outfile)) {
        stop("outfile must be a new intensity backing-file path.")
    }
    if (!dir.exists(dirname(outfile))) stop("outfile must have an existing parent directory.")
    if (file.exists(outfile)) stop("Normalization output already exists; choose a new outfile.")
    n <- .msiNumberOfPixels(x)
    arrays <- methods::is(x, "MSImagingArrays")
    complete <- FALSE
    on.exit(if (!complete && file.exists(outfile)) unlink(outfile), add = TRUE)
    destination <- if (arrays) {
        matter::matter_list(data = NULL, path = outfile, type = "double", readonly = FALSE,
                            lengths = lengths(Cardinal::mz(x)), names = names(Cardinal::intensity(x)))
    } else {
        matter::matter_mat(data = NULL, path = outfile, type = "double", readonly = FALSE,
                           nrow = length(Cardinal::mz(x)), ncol = n,
                           dimnames = dimnames(Cardinal::intensity(x)))
    }
    if (verbose) {
        progress <- utils::txtProgressBar(min = 0, max = n, style = 3)
        on.exit(close(progress), add = TRUE)
    }
    for (pixel in seq_len(n)) {
        spectrum <- .msiSpectrumVectors(x, pixel)
        value <- .msiNormalizeSpectrum(spectrum$intensity, spectrum$mz, method,
            reference, mz, tol, zeroAsMissing, scale)
        if (arrays) destination[[pixel]] <- value else destination[, pixel] <- value
        if (verbose) utils::setTxtProgressBar(progress, pixel)
    }
    result <- x
    Cardinal::intensity(result) <- destination
    Cardinal::pixelData(result)$moleculaRNormalizationFactor <- factors[1, ]
    result <- .msiRecordProcessing(result, "moleculaR intensity normalization",
        list(method = method, scale = scale, meanMethod = meanMethod,
             zeroAsMissing = zeroAsMissing, mz = mz, tol = tol,
             implementation = "Cardinal/matter rescaling; streamed native intensity backing"))
    complete <- TRUE
    result
}

.filterCardinalPeaks <- function(x, minFreq, featureChunkSize) {
    .msiCheckProcessing(x)
    if (!methods::is(x, "MSImagingExperiment")) {
        stop("Peak-frequency filtering requires an aligned MSImagingExperiment; align arrays explicitly first.")
    }
    if (!isTRUE(ProtGenerics::centroided(x))) stop("Peak-frequency filtering requires known centroided features.")
    if (!is.numeric(minFreq) || length(minFreq) != 1L || !is.finite(minFreq) || minFreq < 0 || minFreq > 1) {
        stop("minFreq must be a finite number in [0, 1].")
    }
    count <- numeric(length(Cardinal::mz(x)))
    for (rows in .msiFeatureBlocks(x, featureChunkSize)) {
        values <- as.matrix(Cardinal::intensity(x)[rows, , drop = FALSE])
        .msiValidateIntensities(values)
        count[rows] <- rowSums(values > 0)
    }
    frequency <- count / .msiNumberOfPixels(x)
    Cardinal::featureData(x)$count <- count
    Cardinal::featureData(x)$freq <- frequency
    keep <- which(frequency >= minFreq)
    result <- x[keep, , drop = FALSE]
    .msiRecordProcessing(result, "moleculaR peak frequency filtering",
        list(minFreq = minFreq, occurrence = "positive intensity per shared feature row",
             featuresBefore = length(count), featuresAfter = length(keep)))
}
