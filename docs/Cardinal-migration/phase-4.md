# Phase 4: native preprocessing and measured FWHM

The native scientific workflow now uses Cardinal/matter for spectral
preprocessing and peak picking. Numeric profiles, native precomputed widths and
Cardinal profile objects feed the existing half-height estimator and curve
fits. Default instrument calibrations and the spatial probability algorithms
are preserved. Legacy input adapters remain temporarily available for existing
Shiny/example workflows; removing those dependencies is a later phase.

The user explicitly authorized reasonable Cardinal alternatives for
normalization and picking. Native defaults therefore follow Cardinal's scaling
and zero-reference behavior, rather than enforcing historical implementation
details. Differences below are measured and intentional.

## Input provenance

**All real peak-picking and measured-width evaluation uses only the separate
`inst/extdata/pos-XIII-82492-singleSpectrum.tsv` profile.** Its MD5 is
`73560b65278d3431750fa0884a12b76b`. The supplied imzML is reduced/shrunk
centroided data and is never used to pick peaks or estimate widths.

The centroided imzML is used only after width fitting, for extraction/map
sensitivity and file-backed normalization checks. Those results apply to its
provided reduced feature set. Generated profile imzML fixtures test the native
profile adapter independently of that acquisition. Known-centroided and
unknown-status Cardinal inputs are rejected for measured-width estimation.

## Use

```r
# A separate measured profile; specify whether the TSV has a header.
profile <- readSingleSpect("profile.tsv", header = FALSE)
set.seed(4)
measured <- estimateFwhm(profile, peakMethod = "diff", SNR = 5)

# Equivalent numeric input, or externally supplied measured widths.
# measured <- estimateFwhm(mz = mass, intensity = signal)
# measured <- estimateFwhm(data.frame(mz = peakMasses, fwhm = peakWidths))

# Centroided MSI is suitable for normalization and extraction, not measured FWHM.
msData <- readMsi("dataset.imzML", memory = FALSE)
msData <- normalizeIntensity(msData, "TIC", outfile = "normalized-intensity.bin")
analytes <- batchAnalyteSearch(msData, measured, m = c(760.5851, 788.5447))
if (analytes$n > 0) {
    mpm <- probMap(analytes, seed = 1)
}

# With no separate profile, the unchanged approximate instrument model works.
approximate <- useDefaultFwhm(modality = "FTICR")
```

`readSingleSpect()` returns a numeric `moleculaRSpectrum` list with `mz` and
`intensity`. The default `header = TRUE` preserves the old first-row handling;
use FALSE for headerless files, including the bundled profile. The baseline
comparison retains header=TRUE for both readers, ensuring identical inputs.
Exactly two columns, positive strictly increasing masses, and finite
non-negative intensities are required. Empty, duplicate/unsorted, malformed and
non-finite profiles fail explicitly.

For a Cardinal **profile** object, select one run and complete queued processing
before `estimateFwhm()`. Only sampled spectra are retrieved. `spectraSampling`
and `peakSampling` control sampling through R's random stream; set.seed()
provides reproducibility. Use `peakIndices` for one numeric profile to isolate
width calculation from automatic selection. Models record selected spectra,
picker/settings and detected/discarded counts in `estimation`. Unbounded
calculated widths are discarded with a warning; at least five distinct valid
width samples are required. Failed estimation never substitutes an instrument
model silently.

## Default calibration and measured-width evidence

`inst/extdata/fticr-fwhm-reference.rds` is a 20,564-byte native numerical
calibration: interpolation knots/options, measured values and source provenance.
It contains no function environment or MALDIquant spectrum. `useDefaultFwhm()`
reconstructs the original interpolator without loading the legacy MSI dataset.
Predictions agree with the frozen baseline at knots, a dense grid, custom
resolving-power anchors and extrapolation points. The observed maximum default
prediction difference is zero. The original bundled dataset is retained.

The numeric half-height implementation uses the original two linear
interpolations. On the controlled 21-peak Gaussian profile, supplied-index widths
and the fitted curve reproduce legacy values exactly. Automatic Cardinal diff
picking recovers those same 21 peaks. Native precomputed tables and both
Cardinal profile representations are covered, including backed profile input.

On the separate real TSV profile, at SNR 5:

| Measurement | Result |
|---|---:|
| Legacy SuperSmoother detections | 3,718 |
| Cardinal diff detections | 1,462 |
| Shared profile-grid peak positions | 933 |
| Legacy-only positions | 2,785 |
| Cardinal-only positions | 529 |
| Maximum width difference at common positions | 0 |
| Widths retained for the seeded native fit | 1,000 |
| Invalid/unbounded native widths discarded | 0 |

The automatic fitted curves differ because their peak sets and seeded samples
differ. This is not a half-height calculation discrepancy:

| m/z | Legacy fitted width | Cardinal fitted width | Relative change |
|---:|---:|---:|---:|
| 150 | 0.0006603321 | 0.0006428217 | -2.65% |
| 249 | 0.0018441464 | 0.0020418448 | +10.72% |
| 400 | 0.0046845895 | 0.0040233403 | -14.12% |
| 788.5447 | 0.0162770751 | 0.0162940997 | +0.10% |
| 1,000 | 0.0244392557 | 0.0275009741 | +12.53% |

These are five sampled predictions, not global error bounds or evidence that
either automatic curve is an instrument calibration standard. The default
FTICR instrument model is independent of this comparison and remains unchanged.

Applying the two fitted widths at m/z 788.5447 to the **provided reduced
centroided MSI** gives 17,642 detections for both. Gaussian intensities change in
17,642 pixels, with maximum absolute change 6,815.834 in the supplied intensity
units. With the same tissue domain, bandwidth 2 and seed 1, maximum analyte-density
change is about 1.248146e-7 and three hotspot/coldspot/neutral grid classifications
change. This checks sensitivity on the available reduced data, not a complete
original spectrum or broad biological-impact assessment. No normalization is
mixed into this width-only comparison.

Exact evidence is in `phase4-peak-comparison.csv`, `phase4-fwhm-comparison.csv`
and `phase4-map-comparison.csv`. Reproduce with `validate-phase4.R`.

## Normalization

Native TIC, RMS and reference workers call the exported matter functions used
by Cardinal. Their results are compared directly with Cardinal, while a small
guard handles empty spectra that otherwise fail Cardinal RMS. Queued processing
must be completed first; inputs and nonzero-spectrum factors must be finite.
Empty/all-zero acquired spectra stay zero. Repeated TIC/RMS calls are stable.

| Method | Native behavior |
|---|---|
| TIC | Cardinal sum rescaling; default target sum is the maximum spectrum length/shared feature count |
| RMS | Cardinal RMS rescaling with the same default scale; use scale=1 for unit RMS |
| IS | Cardinal reference rescaling within an explicit absolute tol; absent/zero standards yield zero spectra |
| median | Custom division by the per-spectrum median; undefined nonzero-spectrum factors fail |
| medianFoldChange | Custom aligned feature reference using mean/median/sum, followed by median ratio division |

`scale = 1` reproduces division-only TIC/RMS on equivalent numerical spectra.
Default Cardinal scaling multiplies these results by a dataset-wide constant;
a seeded map test verifies that this constant cancels in normalized density
calculations. It still changes reported intensity values. IS deliberately omits
legacy median imputation, so missing standards can remove analyte detections.
The absolute tolerance prevents normalization against an unrelated distant mass.
Reference normalization requires sorted m/z values and follows Cardinal's tie
behavior.

Native medianFoldChange requires an explicitly aligned, known-centroided
`MSImagingExperiment` and one run. With `zeroAsMissing = TRUE`, shared zeros
represent missing centroid peaks; with FALSE, every shared feature is treated
as observed. Both definitions are tested against matching legacy fixtures.
Reference-zero features are excluded from ratios. Arrays are never aligned
implicitly. Feature-row workspaces use Cardinal's byte setting or a 16 MiB
target, with `featureChunkSize` as an explicit override.

Normalization returns the same native Cardinal class with a new double-precision
intensity backing file. Original m/z storage, coordinates and metadata are kept.
Factors are recorded in `pixelData()$moleculaRNormalizationFactor` and completed
history in processingData. A missing reference factor can be NA in provenance;
output intensities remain finite. Source files are not rewritten. Keep the new
file **and original m/z sources** accessible when persisting the object, or
export it with Cardinal's writer. Existing output paths are refused.

The file write pass holds one retrieved spectrum at a time. This avoids an
installed Cardinal/matter failure encountered while recursively combining 2,000
file-backed process() results. It does not modify Cardinal or use legacy
normalization math. BiocParallel can handle the scalar factor pass; writing is
serial. Worker startup temporarily exposes the isolated library through R_LIBS,
then restores the caller's environment; caller-started backends remain running.

## Filtering

`filterPeaks()` requires aligned centroided native data and a completed queue.
It recomputes positive-intensity pixel counts/frequencies, retains
`freq >= minFreq` including equality, and records counts, frequencies and
processing history. Duplicate shared rows remain distinct. No acquired pixel
is removed, and an empty feature result remains a valid native object.
The same bounded feature-block workspace is used.

Stored zeros count as absent here. Legacy sparse stored-entry counts could
count zero-valued stored peaks, which is a representation difference rather
than an inferred alignment rule. Native arrays require explicit alignment.

## Validation, performance and remaining work

The source suite passes **666 assertions across 40 tests**, with zero failures,
errors, skips or test warnings. The final source archive builds and installs;
R CMD check also passes all 666 installed-package assertions and all four older
standalone regression scripts, with **zero errors, 10 existing warnings and
2 notes**. Remaining categories concern legacy filenames/encoding/imports, S3
signatures, older documentation, data compression and vignette packaging.
PDF/manual generation and vignette rebuilding were disabled; suggested
rmarkdown was unavailable. Full output is in package-build.txt/package-check.txt.
Tests cover numeric import, width interpolation and fitting, precomputed widths,
profile status/queues/run selection, direct Cardinal normalization equivalence,
legacy division-only/custom-reference comparisons, zeros/invalid data, file
immutability, cache invalidation, inclusive frequencies, duplicate features,
empty results, backend ownership and the native path through seeded maps.
Typed Cardinal byte budgets and ordinary record counts have separate checks;
an explicit pixel/feature chunk size retains its stated units.

A local 2,000-spectrum subset of the supplied reduced centroided data normalized
with native TIC in 4.24 s. Output was file-backed: 34,882,968 new binary bytes
and a 626,920-byte native object. Rprofmem recorded no individual R allocation
of at least 1 MiB during normalization. Whole-process peak working set was
1,546,235,904 bytes; this includes real-profile analysis, import and both map
comparisons and is not incremental normalization RAM. Native allocations and
total allocation volume are not measured. This is one local observation with
uncontrolled disk cache/background load, not a large-scale acceptance result.
TIC in this storage check sums only the available reduced features; it does
not reconstruct the original acquisition's full-spectrum total ion current.
Raw figures are in `phase4-benchmark.csv`.

Phase 5 remains Cardinal integration into Shiny. Example/vignette cleanup,
legacy dependency removal, broader biological data, other platforms and large
benchmarks remain later phases. No shrunk centroided acquisition is treated as
a profile, and no probability-map algorithm is changed in this phase.
