# Scientific and architectural decisions

## 2026-11-23 — D001: migration boundary and temporary legacy support

Cardinal owns spectral data; existing moleculaR point patterns and `probMap()`
own spatial analysis. Add a Cardinal path to `searchAnalyte()` while preserving
the legacy sparse path temporarily for regression and existing callers. Do not
change `probMap()`, `gaussWeight()`, the FWHM models, or preprocessing in this
milestone. Remove legacy dependencies only in the validated cleanup phase.

## 2026-11-23 — D002: mass selection

For Cardinal input select measured masses with `abs(mz - m) <= 3 * sigma`,
where `sigma = fwhm / (2 * sqrt(2 * log(2)))`. Both endpoints are inclusive.
Do not add peaks outside this interval by rounding to nearest global indices.
Legacy nearest-index selection is retained only on the temporary sparse path.
Boundary fixtures must quantify the resulting difference. No extra tolerance
or normalization is applied. Gaussian weights use existing `gaussWeight()`.

## 2026-11-23 — D003: explicit mean semantics

For `MSImagingArrays`, mean uses the number of measured peaks in that spectrum's
selected interval (including explicitly stored zeros). An empty interval returns
zero. For `MSImagingExperiment`, mean uses the number of selected shared features,
including zeros, preserving its aligned-data interpretation. This may differ
between representations when absent peaks become zero-filled shared features.
Reconstructing a global feature axis for unaligned means would defeat the target
architecture. Document and test this intentional scientific/API difference.

## 2026-11-23 — D004: processing and input validity

Require callers to execute queued processing with `Cardinal::process()` before
extraction. `spectrapply()` otherwise ignores those operations in Cardinal 3.14.0.
Inspect the queue using `ProtGenerics::processingData()` (not exported by
Cardinal). The centroid-status accessor is also `ProtGenerics::centroided()`.
Never process implicitly, reset queues, or peak-pick centroided data.
Require known centroided status for extraction; import also accepts profile or
unknown status for explicit subsequent processing. Reject non-finite or negative
intensities and non-finite/non-positive masses rather than silently imputing or
discarding them. Permit duplicate or unsorted measured masses; aggregate every
selected observation. Reject invalid/ambiguous spatial coordinates, multiple
runs, and varying spatial z for a single two-dimensional point pattern.

## 2026-11-23 — D005: spatial domain and annotations

Keep spectrum order and coordinates together. Construct a common tissue window
from all acquired pixels, including analyte zeros. Use existing spatial window
construction without physical scaling. Users can supply a valid window when
coordinate spacing cannot be inferred; validate that it contains all acquired
pixel centres. Omit zero-intensity analyte points and preserve the legacy empty
`ppp` return convention. Annotation matching uses nearest absolute mass within
the same tolerance, breaking ties by original input order; `NA` scalar means
no supplied annotations and an empty numeric vector means no confirmations.

## 2026-11-23 — D006: environment and dependency preference

Use the stable Bioconductor 3.23 / R 4.6 environment and Cardinal 3.14.0.
Keep installation inside `renv/library/cardinal-migration` and leave user/global
libraries and the old lockfile untouched until the environment is validated.
After `BiocManager::valid()` returned TRUE, updated the root lockfile and pinned
Bioconductor source records to the official 3.23 release. The pre-migration
lockfile remains available at the preserved starting commit.
Prefer Bioconductor packages where suitable alternatives exist; use BiocParallel
serial execution initially, and ProtGenerics for its exported processing accessor.
No direct matter dependency is needed for implementation; it is used in tests
and performance verification of file-backed storage. Reference documentation:
[Bioconductor installation](https://bioconductor.org/install/) and
[Cardinal 3.23 package page](https://bioconductor.org/packages/3.23/bioc/html/Cardinal.html).

## 2026-11-23 — D007: observed Cardinal runtime differences

In installed Cardinal 3.14.0, `subsetPixels()` calls two-dimensional indexing
even for arrays, whose subset method requires its first index. Use native
class-appropriate indexing in `.msiSelectPixels()` and expose exact run-name
selection through `readMsi(run = ...)`. Existing arrays can be subset with
`x[Cardinal::run(x) == "sample"]`; experiments use the pixel column index.

Cardinal's default TIC normalization rescales the example spectra, producing
selected sums 3 and 0.75 rather than division-only sums 1 and 0.25. The first
milestone honors already-processed intensities without changing this behavior.
Normalization equivalence remains a separate phase 4 task.

Single-run imzML export/import can replace its run label with the filename.
Test reader metadata against direct Cardinal reading of the same file, rather
than assuming that Cardinal's writer preserves an original in-memory run label.
`readMsi(check = TRUE)` promotes reader warnings to errors: Cardinal otherwise
only warns about corrupt binaries. Default `check = FALSE` keeps Cardinal's
fast import behavior and delays binary read errors until spectral access.

## 2026-11-23 — D008: round-trip fixture precision

The installed CardinalIO writer defaults to `intensity.type = "float32"` and
`mz.type = "float64"`. Writing the double-precision legacy fixture with those
defaults changes spectral inputs and produces small density differences around
1e-10 after an imzML round trip. This is a storage-precision difference, not
an extraction discrepancy. Generate controlled regression pairs explicitly
with `intensity.type = "float64"`; verify imported intensities and retain the
original 1e-10 scientific comparison tolerance. Production import preserves
the input file's precision; it does not rewrite or modify the source data.

## 2026-11-23 — D009: explicit spatstat dependencies and readable source

Retire the release-era dynamic importer. Resolve each used spatial function
against the installed child-package exports and call its owner explicitly:
`spatstat.geom`, `spatstat.explore`, `spatstat.random`, or `spatstat.utils`.
Declare those direct dependencies and remove `import` and the `spatstat`
meta-package from runtime dependencies. Standard generics such as
`summary(window)` retain ordinary S3 dispatch; a direct method call uses its
owner, for example `spatstat.geom::summary.owin(window)`.

This intermediary changes name resolution and source readability only. Preserve
spatial algorithms, parameters, data structures, and saved numerical fixtures.
Verify the qualification pass by comparing parsed expressions after removing
the namespace prefixes, then rerun numerical and spatial regressions. Replay
the immutable legacy source in a test-only environment with its three historical
geometry bindings; do not reintroduce runtime dynamic imports.

Expand comments, argument layouts, and validation blocks in the new Cardinal
helpers and development scripts so their purpose, data flow, memory behavior,
and scientific conventions can be followed without reconstructing the migration.

## 2026-11-25 — D010: batch extraction and output semantics

Phase 3 visits each spectrum once per batch, evaluates all unique (mass, FWHM)
queries there, and stores positive results as a sparse pixel-by-query matrix.
Pixel chunks bound transient spectral/results storage; output memory still
grows with the number of detections. Obtain the default pixel chunk size from
Cardinal's configured chunk size or number of chunks; expose an explicit
`chunkSize` override for reproducible tuning. Never build a global spectral
matrix. Shared experiments retain the union of queried feature rows, including
zeros for their mean denominator. Arrays retain their observed-peak denominator.

Batch results must match repeated single-analyte searches, including inclusive
boundaries and windows overlapping the edge of the measured mass range. Do not
discard a query solely because its centre lies outside that range. Missing
(`NA`) targets are skipped; non-positive/non-finite supplied targets or invalid
model widths are errors. Empty/no-detection batches return an empty ppp. Retain
input order and duplicate query metadata, while deduplicating spectral work.
Lipid metadata follow database row order and the established M-H, M+H, M+Na,
M+K adduct order. Use the same three-sigma confirmation rule as single search;
the legacy FWHM-only database prefilter must not discard valid confirmations or
duplicate database rows. The temporary legacy APIs remain unchanged.

Serial execution remains the default. Expose BiocParallel configuration and
map explicit `numCores > 1` to a cross-platform SnowParam, with equivalence tests
and timing comparisons before recommending it. Caller-owned backends retain
their lifecycle; stop only backends started internally for this call.

## 2026-11-25 — D011: optional bounded extraction cache

An optional caller-owned environment retains only the last extraction and its
key. The key includes the native object (spectra, coordinates, run/ROI,
representation, metadata and processing state), requested masses and evaluated
FWHM values, aggregation and mass-selection rule, and backing-file paths,
sizes and modification times. Query annotations and lipid labels are reapplied
to cached intensities. Do not cache generated point-pattern IDs or spatial maps.
Use matter's exported path accessor for backing-file provenance and declare
matter as a direct dependency. Check file state again before saving a result.

Backing data must remain immutable while in use. File statistics do not detect
external rewrites that preserve both size and timestamps; callers must clear or
replace the cache after such edits. No full-file hashing or spectral conversion
is performed to create a cache key. Changing normalization, run/ROI, target,
width, representation, ordinary file timestamps or intensity values invalidates
the entry. The cache is opt-in, session-local and not a durable artifact.
