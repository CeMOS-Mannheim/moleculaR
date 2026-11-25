# Cardinal batch extraction. Spectra remain native and file-backed. Only sparse
# query results are retained after each pixel chunk has been read once.

.msiBatchBackend <- function(numCores, BPPARAM) {
    if (!is.numeric(numCores) || length(numCores) != 1L || !is.finite(numCores) ||
        numCores < 1 || numCores != floor(numCores)) {
        stop("numCores must be a positive integer.")
    }
    if (!is.null(BPPARAM)) {
        if (!methods::is(BPPARAM, "BiocParallelParam")) {
            stop("BPPARAM must be a BiocParallelParam object.")
        }
        if (numCores != 1L) {
            stop("Specify either numCores > 1 or BPPARAM, not both.")
        }
        return(BPPARAM)
    }
    if (numCores == 1L) BiocParallel::SerialParam() else
        BiocParallel::SnowParam(workers = numCores, progressbar = FALSE)
}

.msiBatchChunkSize <- function(n, chunkSize, x = NULL) {
    if (is.null(chunkSize)) {
        chunkSize <- Cardinal::getCardinalChunksize()
        if (length(chunkSize) != 1L || is.na(chunkSize)) {
            chunkSize <- ceiling(n / Cardinal::getCardinalNChunks())
        } else if (inherits(chunkSize, "size_bytes")) {
            if (is.null(x)) stop("A Cardinal object is required for a byte-sized chunk budget.")
            peaks <- if (methods::is(x, "MSImagingArrays")) max(lengths(Cardinal::mz(x))) else length(Cardinal::mz(x))
            # Cardinal accepts both record counts and typed byte budgets. This
            # conservative estimate covers double mass/intensity vectors for
            # the longest spectrum; output memory still grows with detections.
            chunkSize <- max(1, floor(as.numeric(chunkSize) / (16 * max(1, peaks))))
        }
    }
    if (!is.numeric(chunkSize) || length(chunkSize) != 1L || !is.finite(chunkSize) ||
        chunkSize < 1 || chunkSize != floor(chunkSize) || chunkSize > .Machine$integer.max) {
        stop("chunkSize must be a positive integer number of pixels.")
    }
    min(n, chunkSize)
}

# Observe backing files without reading the spectra or hashing complete binaries.
# In-memory objects have no files. Cached data require immutable backing files.
.msiBackingState <- function(x) {
    arrays <- ProtGenerics::spectraData(x)
    paths <- unique(unlist(lapply(seq_along(arrays), function(i) {
        a <- arrays[[i]]
        if (methods::is(a, "matter")) as.character(matter::path(a)) else character()
    }), use.names = FALSE))
    paths <- sort(normalizePath(paths, winslash = "/", mustWork = TRUE))
    info <- file.info(paths)
    data.frame(path = paths, size = info$size, mtime = as.numeric(info$mtime))
}

# Precompute sorted unique queries and the map back to the caller's order.
.msiBatchQueries <- function(mz, fwhm, method) {
    if (!is.numeric(mz) || !is.numeric(fwhm) || length(mz) != length(fwhm)) {
        stop("mz and fwhm must be numeric vectors of equal length.")
    }
    .msiValidateTarget(1, 1, method)
    tolerance <- vapply(seq_along(mz), function(i) {
        .msiValidateTarget(mz[i], fwhm[i], method)
    }, numeric(1))
    ordering <- order(mz, fwhm)
    pairs <- data.frame(mz = mz[ordering], fwhm = fwhm[ordering])
    first <- !duplicated(pairs)
    map <- integer(length(mz))
    map[ordering] <- cumsum(first)
    selected <- ordering[first]
    mz <- mz[selected]
    tolerance <- tolerance[selected]

    # Conservative lookup bounds avoid losing an endpoint through subtraction
    # rounding. The worker always reapplies the exact single-search predicate.
    guard <- 2 * .Machine$double.eps * pmax(abs(mz), tolerance)
    list(mz = mz, fwhm = fwhm[selected], tolerance = tolerance,
         lower = mz - tolerance - guard, upper = mz + tolerance + guard, map = map)
}

# A worker returns only positive query values. Sort an unaligned spectrum once,
# then locate each interval by binary search; duplicates and stored zeros remain.
.msiAggregateIons <- function(intensity, mz, queries, method) {
    .msiValidateMasses(mz)
    if (!is.numeric(intensity) || length(intensity) != length(mz) ||
        any(!is.finite(intensity)) || any(intensity < 0)) {
        stop("Retrieved intensities must be finite, non-negative, and match the m/z values.")
    }
    if (is.unsorted(mz)) {
        ordering <- order(mz, method = "radix")
        mz <- mz[ordering]
        intensity <- intensity[ordering]
    }
    left <- findInterval(queries$lower, mz, left.open = TRUE) + 1L
    right <- findInterval(queries$upper, mz)
    values <- numeric(length(queries$mz))
    for (j in which(left <= right)) {
        indices <- seq.int(left[j], right[j])
        selected <- indices[abs(mz[indices] - queries$mz[j]) <= queries$tolerance[j]]
        if (!length(selected)) next
        signal <- intensity[selected]
        values[j] <- switch(method,
            Gaussian = sum(signal * gaussWeight(mz[selected], queries$mz[j], queries$fwhm[j])),
            sum = sum(signal), mean = mean(signal), max = max(signal))
    }
    if (any(!is.finite(values))) stop("Aggregated analyte intensity is not finite.")
    positive <- which(values > 0)
    list(j = positive, x = values[positive])
}

.msiExtractIons <- function(x, mz, fwhm, method = "Gaussian", chunkSize = NULL,
                           BPPARAM = BiocParallel::SerialParam(), cache = NULL,
                           verbose = FALSE) {
    .msiValidate(x, extraction = TRUE)
    queries <- .msiBatchQueries(mz, fwhm, method)
    n <- .msiNumberOfPixels(x)
    chunkSize <- .msiBatchChunkSize(n, chunkSize, x)
    if (!methods::is(BPPARAM, "BiocParallelParam")) {
        stop("BPPARAM must be a BiocParallelParam object.")
    }
    if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
        stop("verbose must be TRUE or FALSE.")
    }
    if (!is.null(cache) && (!is.environment(cache) || environmentIsLocked(cache))) {
        stop("cache must be NULL or an unlocked environment.")
    }

    # Store one request only. Native object identity includes in-memory values,
    # processing metadata and run/ROI order; file state covers backed spectra.
    key <- if (!is.null(cache)) list(object = x, files = .msiBackingState(x),
        mz = mz, fwhm = fwhm, method = method, rule = "inclusive-three-sigma-v1")
    if (!is.null(cache) && exists(".moleculaR_ions", cache, inherits = FALSE)) {
        previous <- get(".moleculaR_ions", cache, inherits = FALSE)
        if (is.list(previous) && identical(previous$key, key)) return(previous$values)
    }
    result <- Matrix::sparseMatrix(i = integer(), j = integer(), x = numeric(),
                                   dims = c(n, length(mz)))
    if (!length(mz)) return(result)

    if (methods::is(x, "MSImagingExperiment")) {
        axis <- Cardinal::mz(x)
        .msiValidateMasses(axis)
        selected <- rep(FALSE, length(axis))
        for (j in seq_along(queries$mz)) {
            selected <- selected | abs(axis - queries$mz[j]) <= queries$tolerance[j]
        }
        if (!any(selected)) return(result)
        x <- x[which(selected), , drop = FALSE]
    }

    # Reuse a caller-started backend. Stop only a cluster that this call starts.
    started <- !BiocParallel::bpisup(BPPARAM)
    if (started) {
        .msiStartBackend(BPPARAM)
        on.exit(BiocParallel::bpstop(BPPARAM), add = TRUE)
    }
    starts <- seq.int(1, n, by = chunkSize)
    rows <- columns <- signals <- vector("list", length(starts))
    if (verbose) {
        progress <- utils::txtProgressBar(min = 0, max = n, style = 3)
        on.exit(close(progress), add = TRUE)
    }
    for (k in seq_along(starts)) {
        pixels <- seq.int(starts[k], min(n, starts[k] + chunkSize - 1))
        block <- .msiSelectPixels(x, pixels)
        values <- Cardinal::spectrapply(block, FUN = .msiAggregateIons,
            index = "mz", queries = queries, method = method, simplify = FALSE,
            verbose = FALSE, BPPARAM = BPPARAM)
        if (length(values) != length(pixels)) {
            stop("Cardinal batch extraction did not return one result per pixel.")
        }
        counts <- vapply(values, function(value) length(value$j), integer(1))
        rows[[k]] <- rep.int(pixels, counts)
        columns[[k]] <- unlist(lapply(values, `[[`, "j"), use.names = FALSE)
        signals[[k]] <- unlist(lapply(values, `[[`, "x"), use.names = FALSE)
        if (verbose) utils::setTxtProgressBar(progress, tail(pixels, 1L))
    }
    result <- Matrix::sparseMatrix(i = unlist(rows, use.names = FALSE),
        j = unlist(columns, use.names = FALSE), x = unlist(signals, use.names = FALSE),
        dims = c(n, length(queries$mz)))
    result <- result[, queries$map, drop = FALSE]
    if (!is.null(cache)) {
        if (!identical(key$files, .msiBackingState(key$object))) {
            stop("Backing files changed during extraction; the result was not cached.")
        }
        assign(".moleculaR_ions", list(key = key, values = result), envir = cache)
    }
    result
}

# Shared spatial preparation occurs once, before dropping acquired zero pixels.
.msiBatchDomain <- function(x, spwin) {
    .msiValidate(x, extraction = TRUE)
    coordinates <- .msiCoordinates(x)
    if (is.atomic(spwin) && length(spwin) == 1L && is.na(spwin)) {
        spwin <- createSpatialWindow(coordinates, plot = FALSE)
    }
    if (!inherits(spwin, "owin")) stop("spwin must be a spatstat owin window.")
    if (!all(spatstat.geom::inside.owin(coordinates$x, coordinates$y, spwin))) {
        stop("spwin must contain every acquired pixel centre.")
    }
    list(coordinates = coordinates, window = spatstat.geom::as.polygonal(spwin))
}

# Missing targets are placeholders, not fabricated masses. Model failures for
# a real target are errors even when that target eventually has no signal.
.msiBatchTargets <- function(m, fwhmObj, method, verifiedMasses, confirmedOnly) {
    if (!is.numeric(m) || any(is.nan(m)) || any(!is.na(m) & (!is.finite(m) | m <= 0))) {
        stop("m must be numeric, with positive finite masses or NA placeholders.")
    }
    if (!is.logical(confirmedOnly) || length(confirmedOnly) != 1L || is.na(confirmedOnly)) {
        stop("confirmedOnly must be TRUE or FALSE.")
    }
    .msiNearestMass(1, verifiedMasses, 0) # Validate annotations even for empty input.
    indices <- which(!is.na(m))
    if (!inherits(fwhmObj, "fwhm")) stop("fwhmObj must be an S3 object of type 'fwhm'.")
    widths <- if (length(indices)) getFwhm(fwhmObj, m[indices]) else numeric()
    .msiBatchQueries(m[indices], widths, method)
    confirmed <- vapply(seq_along(indices), function(i) {
        !is.na(.msiNearestMass(m[indices[i]], verifiedMasses,
                              .msiValidateTarget(m[indices[i]], widths[i], method)))
    }, logical(1))
    keep <- if (confirmedOnly) confirmed else rep(TRUE, length(indices))
    list(indices = indices[keep], mz = m[indices][keep], fwhm = widths[keep],
         confirmed = confirmed[keep])
}

.msiBatchPatterns <- function(values, targets, domain, metadata = NULL) {
    patterns <- lapply(seq_along(targets$mz), function(j) {
        # A sparse column is expanded only for one analyte at a time.
        signal <- as.numeric(values[, j])
        detected <- which(signal > 0)
        if (!length(detected)) return(NULL)
        extra <- if (is.null(metadata)) list() else metadata[[targets$indices[j]]]
        analytePointPattern(x = domain$coordinates$x[detected],
            y = domain$coordinates$y[detected], win = domain$window,
            intensity = signal[detected], mzVals = targets$mz[j],
            metaData = c(list(mzConfirmed = targets$confirmed[j]), extra))
    })
    patterns <- Filter(Negate(is.null), patterns)
    if (!length(patterns)) return(.craeteEmptySpp(domain$window))
    superimposeAnalytes(patterns, spWin = domain$window, check = FALSE)
}

.searchCardinalAnalytes <- function(spData, fwhmObj, spwin, m, numCores, wMethod,
                                    verifiedMasses, confirmedOnly, verbose,
                                    chunkSize, BPPARAM, cache, metadata = NULL) {
    domain <- .msiBatchDomain(spData, spwin)
    targets <- .msiBatchTargets(m, fwhmObj, wMethod, verifiedMasses, confirmedOnly)
    backend <- .msiBatchBackend(numCores, BPPARAM)
    values <- .msiExtractIons(spData, targets$mz, targets$fwhm, wMethod,
        chunkSize, backend, cache, verbose)
    .msiBatchPatterns(values, targets, domain, metadata)
}

# Database processing is independent of spectral access. Preserve row order and
# the historical adduct order, including distinct annotations at repeated masses.
.msiLipidQueries <- function(sldb, adduct) {
    columns <- c("M-H" = "Exact m/z of [M-H]-", "M+H" = "Exact m/z of [M+H]+",
                 "M+Na" = "Exact m/z of [M+Na]+", "M+K" = "Exact m/z of [M+K]+")
    if (!is.character(adduct) || anyNA(adduct)) stop("adduct must be a character vector.")
    adduct[adduct == "M+k"] <- "M+K"
    if (any(!adduct %in% names(columns))) stop("Unsupported adduct.")
    adduct <- names(columns)[names(columns) %in% adduct]
    required <- c("Lipid ID", "Formula (pH7.3)", "Abbreviation*", "numDoubleBond",
                  "lipidGroup", "chainLength", unname(columns[adduct]))
    if (!is.data.frame(sldb) || !all(required %in% names(sldb))) {
        stop("sldb must contain lipid metadata and the requested adduct mass columns; use loadSwissDB().")
    }
    if (!all(vapply(sldb[unname(columns[adduct])], is.numeric, logical(1)))) {
        stop("Requested lipid mass columns must be numeric.")
    }
    rows <- rep(seq_len(nrow(sldb)), each = length(adduct))
    adducts <- rep(adduct, times = nrow(sldb))
    masses <- vapply(seq_along(rows), function(i) sldb[[columns[adducts[i]]]][rows[i]], numeric(1))
    metadata <- lapply(seq_along(rows), function(i) {
        row <- rows[i]
        list(mode = if (adducts[i] == "M-H") "negative" else "positive",
             adduct = adducts[i], lipidID = sldb[["Lipid ID"]][row],
             sumformula = sldb[["Formula (pH7.3)"]][row], abbrev = sldb[["Abbreviation*"]][row],
             numDoubleBonds = sldb$numDoubleBond[row], lipidClass = sldb$lipidGroup[row],
             chainLength = sldb$chainLength[row])
    })
    list(mz = masses, metadata = metadata)
}
