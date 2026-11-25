# Phase 5: Cardinal Shiny integration

Implemented on Cardinal-dev, starting at 8b61395. Decisions D021-D025 are in
scientific-decisions.md; implementation and validation steps are in change-log.md.

| Files | Responsibility |
|---|---|
| R/shiny-cardinal-data.R | Owned files, parsing, .shinyInitialize, native .shinyAnalyze, lipid grouping and .shinySaveAnalysis |
| R/shiny-cardinal.R | Shared .cardinalAnalysisUi/Server, plotting and application factories |
| R/runShiny.R; inst/shiny/package-app.R; inst/shiny/example-app.R | Validated launcher and thin native app entry points |
| inst/extdata/shiny-example; generate-shiny-example.R | Native crop, plain metadata and deliberate generator |
| tests/testthat/test-cardinal-shiny.R; tests/shinyFwhm.R | Native/session regression tests and installed-app smoke tests |
| DESCRIPTION; man/runShiny.Rd; README.md; inst/shiny/www/README.md | Source collation, API help and user guidance |
| validate-phase5.R; phase5-benchmark.csv; test-results.csv; package-build/check.txt | Reproducible validation and measured evidence |

## Native application workflows

Both `runShiny("package-app")` and `runShiny("example-app")` use the same namespaced
analysis module. Main and Example each own a native file-backed MSImagingArrays,
FWHM model, tissue window, input/settings snapshot and last-query extraction cache.
There are no global datasets, legacy sparse copies or shared spectral axes.
The launcher validates its two supported app names.

Upload one matching known-centroided imzML/ibd pair. The module copies the files
into a unique marked session directory and works from those copies. Choose one
run, optional rectangular ROI and pixel spacing when grid inference is ambiguous.
All selected acquired pixels define the domain, including zero-signal pixels and
holes. Profile/unknown-status MSI is rejected with instructions for explicit
Cardinal preprocessing outside the app. No imzML peak picking or width estimation
is implicit in initialization.

FWHM defaults remain FTICR 123000, TIMS-TOF 40000 and TOF-reflector 10000 at m/z
249, with optional custom resolving power/reference mass. The alternative is a
**separate profile TSV**, with explicit header, Cardinal picker and SNR controls.
All valid selected widths are fitted without random peak subsampling. A failed
measured fit reports an error. Reduced centroided MSI cannot supply measured widths.

Optional TIC/RMS/median/internal-standard normalization uses native Phase 4
processing and creates an owned binary intensity file. Every initialization
rereads the original copied source so repeated loading cannot stack normalization.
No normalization is the default; TIC on reduced data describes only its supplied
features. Data/model/preprocessing/annotation changes invalidate readiness, results
and cache; query/statistical changes clear results. Display changes only redraw.

Single and custom-mass collective maps accept exact theoretical masses within the
selected data's measured range. There is no nearest-feature snapping. Lipid class,
(lyso)GPL ion and saturation views select species/adducts before native batch
extraction, preserving identifiers and repeated-mass annotations. Custom mass
lists deduplicate repeated entries. Optional per-analyte standardization uses the
unchanged z-score/range transformation; constant or insufficient signals produce
an explicit error instead of invalid numbers. Confirmation accepts a pre-filtered
mz CSV or a METASPACE export with two metadata rows.

Map calculations call the unchanged probMap implementation with explicit seed 1
and automatic or fixed bandwidth. Computation runs in action handlers with staged
progress and captured warnings. It is synchronous within each R worker; this
stage does not add asynchronous scheduling for simultaneous users. The default
upload limit is 48 GiB when shiny.maxRequestSize is unset; a configured value is
respected. That limit is not a tested dataset size or resource guarantee.

## Example and downloads

The native example contains 1,024 acquired pixels in a 32-by-32 crop of the
preserved reduced centroided glioblastoma data, exported with float64 intensities.
It includes the original annotation metadata, a 546-row annotation-matched lipid
subset and the original 10,319-row SwissLipids species metadata for Main's default
database. These are preserved reference metadata, not a current database download.
The deliberate generator uses legacy objects only during development; app runtime
reads only the native pair and ordinary data frames. A provenance RDS records the
source checksum, commit, bounds and assumptions. General example/vignette migration
and legacy bundle removal remain Phase 6 work.

RDS downloads contain the point pattern, map, aggregated image, FWHM model,
query/settings, verified masses and package versions. CSV joins each detected
point's coordinates/intensity with analyte/lipid/confirmation metadata. PNG exports
the current probability map. No download serializes a live Cardinal object or
private input paths. Original upload paths can expire after initialization without
breaking analysis. Reset/replacement/session closure drops references before
guarded cleanup; only paths inside that module's marked temporary root are removed.

## Validation and remaining work

Validation covers owned paired/auxiliary files, missing/invalid uploads, profile
and unknown-status rejection, run/ROI/zero-pixel selection, default and measured
models, file-backed normalization without stacking, theoretical masses, single and
collective maps, all lipid grouping modes, constant-signal errors, plotting and
portable downloads. Direct native API comparisons and frozen historical map
fields test the unchanged scientific path. Concurrent mock sessions test ownership
and independent initialization, cache/result invalidation and cleanup. Both
installed app entry points generate native example maps, with separate Main/Example
states and preserved packaged inputs. A live launch returns HTTP 200 with both
namespaced tabs and separate-profile controls.

The final source and installed-package suites each pass **829 assertions across 48 tests**, with zero
failures, errors, skips or test warnings. This includes 163 assertions in eight
new tests, including interrupted-write and uploaded-database lifecycle coverage;
all 666 previous assertions still pass. The final source archive builds, installs
and passes all four standalone scripts, with **zero errors, 10 existing warnings
and 2 notes**. Remaining categories concern legacy filenames/encoding/imports,
S3 signatures, older documentation, data compression and vignette packaging.
Reports are saved in package-build.txt/package-check.txt. PDF/manual generation
and vignette rebuilding were disabled; suggested rmarkdown was unavailable.
Validation uses Windows, R 4.6.1, Cardinal 3.14.0 and matter 2.14.0 in the isolated
project library. Build-version startup notices are separate from test warnings.

The native example workflow initialized in 1.17 s, generated a single MPM in
7.72 s, a two-analyte CPPM in 6.57 s and a 32-analyte lipid CPPM in 7.06 s.
The supplied **reduced centroided** 24,504-pixel dataset initialized, including a
private 741,585,326-byte paired-file copy, in 4.91 s; its MPM took 13.25 s and
two-analyte CPPM 11.61 s. Detected points were 874 and 17,646 respectively.
Both paths stayed file-backed, preserved original file checksums and matched
direct API map fields exactly. Neither dataset was used for picking or widths.
The native experiment objects occupied 327,744 and 7,180,312 bytes. Rprofmem
recorded no individual allocation of at least 1 MB during the example single-map
workflow, and a largest allocation of 2,206,448 bytes for the full supplied data.
Whole-process peaks were 1,351,708,672 and 1,562,378,240 bytes, including package
loading and preceding workflows; these are not incremental application RAM.
Native allocations and total allocation volume are not measured. Runs overlapped
other local validation, so timings are observations with uncontrolled disk cache
and resource contention. See validate-phase5.R and phase5-benchmark.csv. These
figures are not large-scale acceptance evidence or measurements of the original
unreduced acquisition.

There is no installed shinytest2/Chromote or supported browser in
this environment; tests use Shiny's testServer, rendered UI, headless graphics and
a live HTTP smoke check. These do not establish browser interaction or cross-platform
acceptance. PDF/manual generation, vignette rebuilding, general legacy API/data
cleanup and larger biological/scale validation remain later work.
