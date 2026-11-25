# Native measured FWHM adapters. Peak selection belongs to Cardinal; the
# half-height interpolation and curve fitting retain moleculaR's definitions.

.msiIsNumericSpectrum <- function(x) {
    is.list(x) && "mz" %in% names(x) && any(c("intensity", "fwhm") %in% names(x))
}

.estimatePeakWidths <- function(mz, intensity, peakIndices) {
    profile <- .numericSpectrum(mz, intensity)
    n <- length(profile$mz)
    if (!is.numeric(peakIndices) || any(!is.finite(peakIndices)) ||
        any(peakIndices != floor(peakIndices)) || any(peakIndices < 1 | peakIndices > n) ||
        anyDuplicated(peakIndices)) stop("peakIndices must contain unique valid sample indices.")
    vapply(peakIndices, function(i) {
        if (i <= 1 || i >= n || intensity[i] <= 0) return(NA_real_)
        half <- intensity[i] / 2
        left <- right <- i
        while (left > 1 && intensity[left] > half) left <- left - 1L
        while (right < n && intensity[right] > half) right <- right + 1L
        if (intensity[left] > half || intensity[right] > half) return(NA_real_)
        xleft <- stats::approx(intensity[left:(left + 1)], mz[left:(left + 1)], xout = half)$y
        xright <- stats::approx(intensity[(right - 1):right], mz[(right - 1):right], xout = half)$y
        abs(xleft - xright)
    }, numeric(1))
}

.msiPickProfile <- function(profile, peakMethod, SNR, backend) {
    # Only this selected profile is materialized. Never peak-pick the complete
    # MSI experiment merely to estimate a calibration from a few spectra.
    input <- Cardinal::MSImagingArrays(
        spectraData = list(mz = list(profile$mz), intensity = list(profile$intensity)),
        pixelData = Cardinal::PositionDataFrame(data.frame(x = 1, y = 1)),
        centroided = FALSE, continuous = FALSE)
    picked <- Cardinal::process(Cardinal::peakPick(input, method = peakMethod, SNR = SNR),
                                BPPARAM = backend, verbose = FALSE)
    peaks <- as.numeric(Cardinal::mz(picked)[[1]])
    if (!length(peaks)) return(integer())
    left <- pmax(1L, pmin(length(profile$mz), findInterval(peaks, profile$mz)))
    right <- pmin(length(profile$mz), left + 1L)
    indices <- ifelse(abs(profile$mz[right] - peaks) < abs(profile$mz[left] - peaks), right, left)
    sort(unique(indices))
}

.estimateNativeFwhm <- function(s, spectraSampling, peakSampling, numCores, smoothingMethod,
                                SNR, plot, savePlot, mz, intensity, peakIndices,
                                peakMethod, BPPARAM) {
    for (name in c("spectraSampling", "peakSampling")) {
        value <- get(name)
        if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
            value < 1 || value != floor(value)) stop(name, " must be a positive integer.")
    }
    if (!is.numeric(SNR) || length(SNR) != 1L || !is.finite(SNR) || SNR < 0) stop("SNR must be finite and non-negative.")
    smoothingMethod <- match.arg(smoothingMethod, c("loess", "superSmoother"))
    peakMethod <- match.arg(peakMethod, c("diff", "sd", "mad", "quantile", "filter", "cwt"))
    if (!is.logical(plot) || length(plot) != 1L || is.na(plot)) stop("plot must be TRUE or FALSE.")
    backend <- .msiBatchBackend(numCores, BPPARAM)

    if (!is.null(mz)) {
        if (!is.null(s)) stop("Supply either s or mz/intensity, not both.")
        s <- .numericSpectrum(mz, intensity)
    } else if (is.numeric(s)) {
        s <- .numericSpectrum(s, intensity)
    } else if (!is.null(intensity)) {
        stop("intensity requires a numeric s or mz vector.")
    }

    imaging <- .msiIsObject(s)
    if (imaging) {
        .msiCheckProcessing(s)
        if (length(unique(.msiRun(s))) != 1L) stop("Select one run before fitting a measured FWHM model.")
        if (!identical(ProtGenerics::centroided(s), FALSE)) {
            stop("Measured FWHM requires known profile spectra; centroided/unknown data need separately supplied precomputed widths.")
        }
        n <- .msiNumberOfPixels(s)
        read <- function(i) .msiSpectrumVectors(s, i)
    } else {
        inputs <- if (.msiIsNumericSpectrum(s)) list(s) else s
        if (!is.list(inputs) || !length(inputs) ||
            !all(vapply(inputs, .msiIsNumericSpectrum, logical(1)))) {
            stop("s must contain native mz/intensity profiles or mz/fwhm tables.")
        }
        n <- length(inputs)
        read <- function(i) inputs[[i]]
    }
    if (!is.null(peakIndices) && n != 1L) stop("peakIndices applies to one profile; use precomputed mz/fwhm tables for multiple spectra.")
    selected <- if (n == 1L) 1L else sample(seq_len(n), min(n, spectraSampling))
    widths <- vector("list", length(selected))
    discarded <- integer(length(selected))
    detected <- integer(length(selected))
    choices <- character(length(selected))
    started <- !BiocParallel::bpisup(backend)
    if (started) {
        .msiStartBackend(backend)
        on.exit(BiocParallel::bpstop(backend), add = TRUE)
    }
    notified <- FALSE
    for (k in seq_along(selected)) {
        input <- read(selected[k])
        if ("fwhm" %in% names(input)) {
            .msiValidateMasses(input$mz)
            if (!is.numeric(input$fwhm) || !length(input$fwhm) || length(input$fwhm) != length(input$mz) ||
                any(!is.finite(input$fwhm)) || any(input$fwhm <= 0)) {
                stop("Precomputed fwhm values must be positive finite numbers matching mz.")
            }
            keep <- rep(TRUE, length(input$mz))
            if (!is.null(input$snr)) {
                if (!is.numeric(input$snr) || length(input$snr) != length(input$mz) ||
                    any(!is.finite(input$snr)) || any(input$snr < 0)) stop("Precomputed snr must be finite non-negative values matching mz.")
                keep <- input$snr >= SNR
            }
            table <- data.frame(peaks = input$mz[keep], fwhmValues = input$fwhm[keep])
            choices[k] <- "precomputed"
        } else {
            profile <- .numericSpectrum(input$mz, input$intensity)
            if (is.null(peakIndices)) {
                if (!notified) {
                    message("Using Cardinal ", peakMethod, " peak detection for measured FWHM (SNR ", SNR, ").")
                    notified <- TRUE
                }
                indices <- .msiPickProfile(profile, peakMethod, SNR, backend)
                choices[k] <- paste0("Cardinal::peakPick/", peakMethod)
            } else {
                indices <- peakIndices
                choices[k] <- "supplied peakIndices"
            }
            table <- data.frame(peaks = profile$mz[indices],
                fwhmValues = .estimatePeakWidths(profile$mz, profile$intensity, indices))
        }
        detected[k] <- nrow(table)
        valid <- is.finite(table$fwhmValues) & table$fwhmValues > 0
        discarded[k] <- sum(!valid)
        table <- table[valid, , drop = FALSE]
        if (nrow(table) > peakSampling) table <- table[sample(seq_len(nrow(table)), peakSampling), , drop = FALSE]
        widths[[k]] <- table[order(table$peaks), , drop = FALSE]
    }
    table <- do.call(rbind, widths)
    table <- table[order(table$peaks), , drop = FALSE]
    if (sum(discarded)) warning("Discarded ", sum(discarded), " unbounded or invalid measured peak widths.")
    if (length(unique(table$peaks)) < 5L) stop("At least five distinct valid peak widths are required for FWHM fitting.")
    interpolator <- switch(smoothingMethod,
        loess = .loess(table), superSmoother = .superSmoother(table$peaks, table$fwhmValues))
    result <- fwhm(interpolator, table$peaks, table$fwhmValues)
    result$estimation <- list(selectedSpectra = selected, detectedPeaks = detected,
        discardedWidths = discarded, peakSelection = choices, SNR = SNR,
        peakSampling = peakSampling, smoothingMethod = smoothingMethod,
        Cardinal = as.character(utils::packageVersion("Cardinal")))
    .plotFwhmResult(result, plot, savePlot)
    result
}
