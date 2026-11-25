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

## 2026-11-25 — D012: independent FTICR calibration

Extract the original interpolator's numeric x/y knots, interpolation method,
endpoint values and options, together with its measured masses/widths, into a
small native RDS under inst/extdata. Reconstruct it using stats::approxfun and
retain the existing positive multiplicative resolving-power calibration. Verify
reference points, a dense grid and extrapolation against the frozen Phase 3
baseline. Default model construction must no longer load the legacy dataset or
any MALDIquant object. Preserve the original dataset for later cleanup.

## 2026-11-25 — D013: native measured spectra and peak detection

readSingleSpect returns a validated numeric moleculaRSpectrum with mz/intensity
vectors. Its default header=TRUE preserves the old reader's first-row handling;
header=FALSE permits explicitly headerless input. Native profile masses must be
positive and strictly increasing; intensities must be finite and non-negative.
estimateFwhm accepts these spectra, numeric vectors, lists/data frames, selected
Cardinal profile spectra, and native precomputed mz/fwhm tables. Existing legacy
input dispatch remains available.

Preserve the two half-height linear interpolations and the existing loess and
SuperSmoother curve fitting. Explicit peakIndices isolate width equivalence.
For automatic native peak detection use Cardinal peakPick(method="diff") by
default, with an explicit peakMethod option. This replaces MALDIquant's
SuperSmoother detector for new native inputs, not its curve-fitting method.
Report that choice and quantify peak-set, width-curve and downstream differences
on synthetic and bundled real profiles. Never assume the pickers are equivalent.
Reject known-centroided Cardinal inputs unless widths are supplied separately;
unknown status is not evidence of profile data. Reject queued processing.
Discard unbounded/invalid calculated widths with an explicit warning, record
counts and selected spectra/peaks, and require at least five distinct valid
width samples for fitting. No fallback to an approximate model after failure.

## 2026-11-25 — D014: native normalization and processing state

Keep the public normalization methods and division-only legacy formulas.
Cardinal's default TIC/RMS scale must not be inherited silently. Use a scalar
factor pass and a custom native addProcessing worker, then Cardinal::process
with outfile to stream double-precision output into a new imzML bundle.
Processing without outfile can materialize the entire result and is unsuitable.
Refuse existing output bundles/intermediate binaries and never overwrite source
files. Return completed native data, recorded factors in pixelData, and Cardinal
processing history. Pre-existing queues must be processed explicitly first.

Empty/all-zero acquired spectra remain zero with factor one. Undefined factors
for nonzero spectra are errors rather than NaN/Inf. RMS uses a scaled equivalent
formula to avoid intermediate overflow. Internal standards use the nearest
measured peak within inclusive absolute tolerance, with deterministic lower-mass
ties; missing standards use the median positive observed standard. A measured
zero standard in a nonzero spectrum is an error. Do not mix runs when estimating
global standard/reference factors. Retain legacy dispatch unchanged.

## 2026-11-25 — D015: aligned median fold-change semantics

Require an explicitly aligned MSImagingExperiment for native medianFoldChange;
do not align arrays silently. Reference mean/median/sum is calculated feature by
feature in bounded blocks. Default zeroAsMissing=TRUE interprets shared-axis zeros
as absent centroid peaks, matching legacy mergeMassPeaks(ignore.na=TRUE) on
positive observed peaks. Set zeroAsMissing=FALSE to reproduce legacy lists that
explicitly store every aligned feature, including zeros. This representation
ambiguity cannot be inferred from intensities alone and must be documented.
Reference-zero features contain no signal and are excluded from ratios; invalid
nonzero-spectrum median factors are errors. Use an explicit featureChunkSize or
Cardinal's configured byte budget, falling back to a documented 16 MiB target
for dense feature-block workspace (at least one feature row).

## 2026-11-25 — D016: feature-frequency filtering

Require aligned, known-centroided data with a completed processing queue.
Recompute count/freq from finite non-negative intensities: one positive value
per acquired pixel contributes to its feature row. Keep freq >= minFreq,
including equality, following the legacy code rather than its old <= prose.
Retain separate shared feature rows even when their masses coincide. Stored
zero values are absent for frequency purposes; legacy sparse stored-entry counts
could include them. Test and document that ambiguity explicitly. Never perform
implicit alignment for arrays. Preserve metadata, empty feature results, and
processing provenance without moving the spatial-analysis boundary.

## 2026-11-25 — D017: user-authorized Cardinal alternatives

The user explicitly clarified that legacy normalization and peak-picking
specifics need not be retained where reasonable Cardinal alternatives exist.
This supersedes D014's division-only scaling, missing-standard median imputation
and zero-standard error requirements for native inputs. D013's Cardinal picker
remains the native default. Legacy input adapters may remain temporarily without
constraining the new native path.

Use Cardinal's actual exported matter rescale_sum/rescale_rms/rescale_ref
implementations inside a Cardinal processing worker. Default scale follows
Cardinal: the maximum spectrum length/shared feature count for TIC/RMS, and one
for reference normalization. Expose scale=1 when division-only output is wanted.
Guard empty spectra (Cardinal RMS currently errors there), reject non-finite
input/factors, and preserve acquired zero pixels. For IS, retain an explicit
absolute matching tolerance as a safeguard against normalizing to an unrelated
peak. Missing or zero references yield zero spectra, following Cardinal's
zero-reference policy, rather than legacy median imputation. Reference axes must
be sorted; ties follow the installed Cardinal implementation. Quantify these
intentional changes against the frozen legacy fixtures and direct Cardinal.

Median and median fold-change have no equivalent in the installed Cardinal
normalization methods and remain custom native processing workers. Their
documented feature/zero handling still follows D015. Preserve half-height width
calculations and default instrument calibrations; changing peak selection is
explicit and tested, not a claim of identical measured-width fits.

## 2026-11-25 — D018: robust native output and isolated worker startup

The initial Cardinal process(outfile=...) prototype worked on tiny inputs but
failed while recursively combining backed results on a 2,000-spectrum real
subset in this installation. process() without outfile returns resident results.
Replace that output mechanism with a double-precision matter_list/matter_mat
destination inside the native Cardinal object. Write one normalized spectrum at
a time, retain original m/z backing, and record completed processing and factors
through exported accessors. This supersedes D014's imzML-bundle output plan;
outfile now names a new intensity backing file, not an imzML export. Existing
paths are refused; partial new files are removed on failure. Persistent objects
require both their new intensity file and original m/z sources to remain present.
The algorithms and scaling remain the approved Cardinal/matter alternatives.

Snow startup also requires the isolated project library to be visible before
workers can load their dependencies. Expose the parent's .libPaths() through
R_LIBS during internally owned startup only, restoring the original environment
afterward. Test caller-owned backend lifecycle and scalar-factor parallel
equivalence; the file write pass remains serial and bounded by one spectrum.

## 2026-11-25 — D019: distinguish profile and reduced centroided evidence

The user clarified that the supplied imzML data are reduced/shrunk centroided
data and unsuitable for picking or estimating widths. Automatic picker and
measured-width evaluation use only the separate bundled
pos-XIII-82492-singleSpectrum.tsv profile, with input checksum recorded.
The real imzML pair is used only after fitting, for extraction/map sensitivity
and file-backed normalization checks. Those results describe the provided
reduced feature set, not an original complete spectral acquisition. Synthetic
profile imzML fixtures test native adapters independently. Reject known-centroided
and unknown-status Cardinal inputs for measured-width estimation; no assumption
of profile data is made from the processed/continuous file-layout label.

## 2026-11-25 — D020: record counts versus typed byte budgets

The installed Cardinal settings distinguish a plain chunk-size number (records)
from a named byte size such as c(MB=2), returned as size_bytes in actual bytes.
Use only a typed byte setting for feature-block workspace, otherwise the 16 MiB
fallback. The batch pixel helper now also distinguishes typed bytes from counts,
estimating pixels conservatively from the longest input spectrum's double
mass/intensity vectors. Explicit pixel chunkSize continues to override settings.
This corrects the old interpretation of a byte budget as a pixel count. Defaults
with NA and ordinary numeric record counts keep their existing behavior.

## 2026-11-25 — D021: one native Shiny session per analysis tab

Replace the two duplicated app servers with a shared module used by upload and
example tabs. Each module owns one native Cardinal MSI object, model, tissue
window, run/ROI selection, processing snapshot and last-request ion cache.
Never assign datasets or analysis results into the global environment. Copy
paired uploads and optional profile/annotation files into unique owned temporary
directories; derived intensity files live there as well. Retire old references
and cache entries before deleting replaced derivatives, and verify ownership
and resolved path boundaries before recursive cleanup. Closing a session removes
only its own temporary files, never packaged examples or original upload paths.

## 2026-11-25 — D022: measured profiles, theoretical masses and input scope

Support only known-centroided two-dimensional MSI in the apps. Profile or
unknown-status MSI is rejected with instructions for explicit external Cardinal
preprocessing; no app peak picking or width estimation is performed on imzML.
Measured FWHM uses a separately uploaded numeric profile TSV (header choice and
Cardinal picker/SNR shown explicitly); otherwise use the preserved instrument
models. Require one run and allow an explicit rectangular ROI. Derive the tissue
window from every acquired pixel in that selection, retaining zeros and holes.
Optional pixel spacing resolves small/irregular grid ambiguity without physical
unit conversion. Theoretical m/z values are never snapped to a measured/shared
axis. Validate centres against the selected data's supported range by streaming
m/z-only summaries for arrays; never construct a global spectral mass axis.

## 2026-11-25 — D023: explicit preprocessing and current results

Offer no normalization by default, plus native TIC/RMS/median/reference methods.
Every initialization rereads the owned input, then applies the selected settings
once; normalization is never stacked accidentally. A changed run, ROI, profile,
model or preprocessing setting marks the session unready, clears results and
invalidates extraction cache until reinitialized. Query/statistical changes
clear results; presentation changes only redraw existing results. Compute maps
in action handlers, never renderPlot. Keep progress stages visible. Calculations
remain synchronous in each R worker; caching and batched extraction avoid
unnecessary rereads but do not promise asynchronous multi-session scheduling.

## 2026-11-25 — D024: native cropped example and lipid grouping

Generate a small 32-by-32 coordinate crop from the preserved reduced centroided
example as a float64 native imzML pair, plus plain annotation/lipid metadata.
The one-time generator can use legacy classes; neither app loads those classes
or the legacy dataset at runtime. Use default FTICR calibration explicitly and
label the crop reduced, not a profile or complete biological acquisition. Retain
SwissLipids identifiers and all supported adduct metadata. Class, ion-milieu and
saturation views select database rows/adducts before native batch extraction;
keep the existing (lyso)GPL class definition. Preserve separate lipid annotations
at repeated masses. Custom mass-list CPPMs deduplicate masses as the old app did.
Optional per-analyte z-score/range standardization calls the existing function;
reject constant/insufficient groups rather than silently changing that formula.

## 2026-11-25 — D025: reproducible maps and self-contained downloads

Keep probMap and its statistical defaults unchanged; expose automatic bandwidth
or a supplied fixed value, with seed 1 as an explicit reproducible UI default.
Display the selected aggregation honestly: a Gaussian-weighted analyte image is
not labeled as an unweighted summed ion image. Downloads contain point values,
lipid/confirmation metadata, map results, model/settings and source filenames;
they do not serialize a live Cardinal object or temporary dataset paths. Export
CSV, RDS and PNG from the current result only. Capture scientific warnings for
display. Record and test intentional model/ROI/preprocessing differences rather
than claiming equivalence to the old app's nearest-feature snapping or globals.

For the separate uploaded profile, fit all valid widths selected by the chosen
Cardinal picker/SNR, without random peak subsampling. This is bounded to one
profile and makes repeated initialization reproducible. The map seed controls
the existing CSR calculation independently. The packaged full SwissLipids species
metadata are also exported to a plain RDS asset, so Main can retain its default
database without loading the legacy spectral bundle.
