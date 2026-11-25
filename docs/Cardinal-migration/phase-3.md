# Phase 3: Cardinal batch analyte and lipid extraction

Both `batchAnalyteSearch()` and `batchLipidSearch()` now accept native
`MSImagingArrays` and `MSImagingExperiment` objects. Spectra remain native and
file-backed. The temporary sparse-data APIs remain available with their existing
behavior; their legacy mass-window/window-diagnostic defects are not repaired
by this stage.

## Use

In the validated local R 4.6.1 / Cardinal 3.14.0 environment:

```r
msData <- readMsi("dataset.imzML", memory = FALSE)
model <- useDefaultFwhm()
targets <- c(760.5851, 788.5447, 810.5290)

analytes <- batchAnalyteSearch(msData, model, m = targets)
mpm <- probMap(analytes, seed = 1)

# Use the SwissLipids table returned by loadSwissDB().
lipids <- batchLipidSearch(msData, model, sldb = sldb,
    adduct = c("M+H", "M+Na"), verifiedMasses = annotations$mz,
    confirmedOnly = TRUE)

# Optional: retain only the last sparse extraction in a session-local cache.
ionCache <- new.env(parent = emptyenv())
analytes <- batchAnalyteSearch(msData, model, m = targets, cache = ionCache,
                             chunkSize = 1000)
```

Select one run/z plane and explicitly finish queued preprocessing with
`Cardinal::process()` before extraction. Profile/unknown centroid status is
rejected. Use `spwin` for grids whose spacing cannot be inferred. The tissue
window always includes all acquired pixels, including zero analyte signal.
Only positive detections become points. Guard against empty results before
calling `probMap()`.

## Architecture and scientific behavior

`R/msi-extractIons.R` validates and deduplicates (mass, evaluated FWHM) queries,
visits each spectrum once, and evaluates all queries in that spectrum. It sorts
unsorted measured masses once and uses binary interval lookups followed by the
exact single-search predicate `abs(mz - target) <= 3 * sigma`. Gaussian weights
continue to use unchanged `gaussWeight()`. Duplicate measured peaks and stored
zeros are retained. Arrays average measured interval peaks; experiments average
the selected shared features, including zeros. No new mean definition is added.

Spectra and results are processed in pixel blocks. The default block size comes
from Cardinal's configured chunk size or number of chunks; `chunkSize` overrides
the pixel count. Experiments select the union of queried feature rows first.
Results are sparse pixels-by-queries, with original query order and duplicate
queries restored. This matrix is an extraction result, not a replacement MSI
backend or global mass-axis conversion.

NA target placeholders are skipped. Other invalid targets, invalid model widths,
non-finite/negative retrieved intensities, and ambiguous spatial input fail
explicitly. A target outside the mass range is still evaluated if its interval
overlaps measured peaks, matching a single search. Empty or wholly undetected
batches return the existing empty `ppp` convention.

Lipid processing generates theoretical masses and metadata before spectral
access. Database-row order and the established M-H, M+H, M+Na, M+K ordering are
preserved. Distinct lipid annotations at repeated masses remain distinct in
metadata even when spectral work is shared. Confirmation uses the single-search
three-sigma rule. The native path omits the legacy database prefilter, whose
one-FWHM tolerance could discard valid confirmations and repeat database rows.
This difference is documented and tested separately from spectral equivalence.

## Cache and parallel execution

The optional environment stores one extraction request. Its key includes the
native data and processing metadata, representation, run/ROI and coordinate
order, target masses, evaluated FWHM values, aggregation, selection rule, and
backing-file paths/sizes/timestamps. Annotations and lipid labels are reapplied
each time; generated point-pattern IDs and probability maps are never cached.
Backing files must stay immutable. Clear/replace the environment after external
edits that preserve file size and modification time; those edits are not
detectable from file statistics. No full spectral hashing is performed.
Replace the cache with `ionCache <- new.env(parent = emptyenv())` to discard
its entry explicitly.

Serial execution is the default. `numCores = 2` requests cross-platform
`BiocParallel::SnowParam` workers. Alternatively pass `BPPARAM` and keep
`numCores = 1`. A backend already started by the caller stays running; a backend
started during the call is stopped on exit. Worker startup and disk contention
can make parallel extraction slower; measured results below determine local
recommendations. No Unix fork behavior is required on Windows.

## Validation and performance

The regression suite compares every aggregation with independent single
extraction on both native representations, multiple pixel chunk sizes and
file-backed input. It checks boundaries, unsorted/duplicate peaks and targets,
empty spectra/batches, NA targets, out-of-range overlaps, annotations, all
adducts, exact lipid metadata, cache invalidation and annotation refresh,
backend ownership, and serial/Snow equivalence. Instrumented workers assert
exactly one visit per acquired spectrum and no visits on a cache hit. Seeded
batch-derived maps are compared with the immutable legacy density and
hotspot/coldspot references at the original 1e-10 tolerance.

The source suite passes **512 assertions across 28 tests**, with zero failures,
errors, test warnings or skips. Test counts are saved in `test-results.csv`.
The source archive builds and installs successfully. R CMD check of that
archive passes all 512 installed-package assertions and all four standalone
regression scripts, with **zero errors, 10 existing warnings and 2 notes**.
The warnings remain in legacy filenames, plotting encoding/imports, S3 method
signatures, older documentation/data compression and vignette packaging;
this stage adds no package-check warning category. PDF generation and vignette
rebuilding were disabled, and suggested rmarkdown was unavailable. Full output
is in `package-build.txt` and `package-check.txt`.

The final benchmark ran in separate fresh R processes without a simultaneous
test run, on Windows with R 4.6.1, Cardinal 3.14.0 and matter 2.14.0. It used the
read-only local `pos-XIII-82492-msi.imzML` pair: 24,504 acquired pixels and a
681,849,892-byte binary. Eight theoretical lipid masses were evaluated using
the unchanged default FTICR FWHM model. Both full-data methods returned 162,416
positive pixel/query pairs in a 1,950,536-byte sparse result.

| Full-data measurement | Repeated single extraction | Batch extraction | Cache hit |
|---|---:|---:|---:|
| Targets | 8 | 8 | 8 |
| Extraction seconds | 36.81 | 8.15 | 0.03 |
| Import seconds | 2.67 | 2.59 | 2.62 |
| Whole-process peak working set, bytes | 724,180,992 | 579,817,472 | 575,758,336 |
| Largest profiled R allocation, bytes | 1,568,304 | 1,299,376 | 0 |

On these observations the batch was approximately 4.5 times faster than eight
separate single searches. Cache-hit timing excludes its initial extraction;
its peak memory includes that warmup. Every batch, cache and tuning result had
**maximum absolute intensity difference zero** from the corresponding repeated
single reference. No scientific tolerance was relaxed.

Pixel-chunk and worker evaluation on the first 2,000 acquired spectra:

| Pixels per chunk | Workers | Extraction seconds |
|---:|---:|---:|
| 100 | 1 | 4.36 |
| 500 | 1 | 1.39 |
| 2,000 | 1 | 0.92 |
| 2,000 | 2 Snow workers | 2.88 |

Small blocks added substantial overhead on this sample. Worker startup/shutdown
is included, so these measurements favor serial execution and larger blocks
for this local case. They do not establish a universal chunk size or parallel
scaling rule. The default continues to honor Cardinal's configuration.

Whole-process peaks include R/package loading, XML parsing and import; they
are not incremental extraction RAM. The parallel figure in the CSV covers the
parent process only and excludes workers. Rprofmem records individual R
allocations, not native allocations or total memory. Disk cache was not flushed,
no physical disk-byte counter was measured, and no hardware/storage inventory
or repeated statistical performance study was performed. Cardinal's existing
scan-direction ontology warning is suppressed only inside the benchmark's
import; production reader behavior is unchanged. This is one moderate-size
dataset, not Phase 7 large-scale validation or a legacy preprocessing benchmark.

Raw metrics are in `batch-benchmark.csv`. Reproduce with
`Rscript --vanilla docs/Cardinal-migration/benchmark-batch.R [imzML]` after
installing the current checkout into the repository-local migration library.
The script checks all results against its independent single reference before
writing the CSV. Larger binary result artifacts remain in ignored local storage.

## Scope and remaining work

Memory grows with detected pixel/query pairs and returned point patterns.
Sparse triplets and the final sparse result coexist temporarily during
assembly; no constant-memory output claim is made. Dense signals across many
thousands of targets can still require substantial RAM. Single-query dense
vectors are created while constructing point patterns, and spatial algorithms
retain their independent memory costs.

Phase 4 remains measured-FWHM input adapters, normalization equivalence and
feature filtering. Later phases cover Cardinal Shiny integration, examples and
vignettes, legacy dependency removal, other operating systems and larger-scale
validation. Existing default FWHM models, probability mapping and historical
fixtures are preserved in this stage.
