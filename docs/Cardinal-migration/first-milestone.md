# First milestone: native Cardinal input to molecular probability maps

The single-analyte workflow now accepts `MSImagingArrays` and
`MSImagingExperiment` directly. Existing molecular probability mapping and FWHM
models are unchanged. The legacy sparse path remains temporarily available for
validation and existing workflows; the overall migration is not complete.

## Run the workflow

The validated environment is R 4.6.1, Bioconductor 3.23, Cardinal 3.14.0,
matter 2.14.0, CardinalIO 1.10.0, and BiocParallel 1.46.0 on Windows 11.
The root `renv.lock` records exact dependency versions. Restore it with renv
in a compatible R environment. For development, the repository-local setup
script installs into the isolated `renv/library/cardinal-migration` library:

```powershell
& 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe' --vanilla docs/Cardinal-migration/setup-environment.R
```

The setup script installs currently available compatible binaries, while the
lockfile is the exact version reference for reproducing this validation. To use
that library interactively from the repository root:

```r
.libPaths(c(file.path(getwd(), "renv/library/cardinal-migration"), .libPaths()))
pkgload::load_all(".")

msData <- readMsi("dataset.imzML", memory = FALSE)
model <- useDefaultFwhm(modality = "FTICR")
target <- 788.5447
analyte <- searchAnalyte(target, getFwhm(model, target), msData)
mpm <- probMap(analyte, seed = 1)
plot(mpm, what = "MPM")
```

Keep the paired imzML/ibd files accessible while the Cardinal object is in use.
No global sparse matrix or mandatory peak alignment is needed. Arrays retain
individual mass axes; experiments select relevant feature rows before chunked
extraction. `readMsi()` imports profile data, but extraction requires known
centroided status and a completed processing queue. Run `Cardinal::process()`
explicitly after any preprocessing. A processed imzML file does not by itself
establish centroid status.

Use `readMsi(..., run = "run-name")` to select an imported run. For existing
arrays use `msData[Cardinal::run(msData) == "run-name"]`; experiments use
`msData[, Cardinal::run(msData) == "run-name", drop = FALSE]`. This avoids an
observed `subsetPixels()` array indexing issue in Cardinal 3.14.0. Spatial
coordinates retain their original units and orientation. A single spatial
z plane is required. Supply `spwin` when the acquisition grid is too small
or irregular for the existing spatial-window inference.

## Numerical interpretation

- Gaussian, sum, mean, and max are supported. The selected interval is
  inclusive: `abs(mz - target) <= 3 * sigma`, with sigma derived from FWHM.
  The unchanged `gaussWeight()` supplies Gaussian weights.
- Strict interval boundaries can exclude a peak admitted by legacy nearest
  feature indices. A controlled boundary example changes sum from 20 to 12
  by excluding one outside-window peak of intensity 8.
- Arrays average measured peaks for `mean`, including stored zeros.
  Experiments average selected shared features, including zero-filled missing
  peaks. A single observed peak of intensity 3 therefore has mean 3 in an
  array, versus 1 on a three-feature shared axis.
- Empty selections return zero. Zero-valued acquired pixels remain part of
  the tissue window; only positive analyte intensities become points. No
  detections return the existing empty `ppp` convention.
- METASPACE confirmation, `confirmedOnly`, analyte metadata, coordinate order,
  and point-mark identifier relationships are preserved. Invalid retrieved
  values, ambiguous coordinates, multiple runs, and queued processing fail
  explicitly rather than silently changing the analysis.

Detailed rationale is in [scientific-decisions.md](scientific-decisions.md).

## Validation and performance

The saved pre-migration function and compact numerical fixtures are in
`tests/fixtures/legacy/`; the starting commit is preserved on
`legacy/MALDIquant-reference`. The automatic suite contains 267 passing
assertions across 14 tests. Equivalent numerical inputs reproduce legacy
intensities, coordinates, seeded MPM densities, and hotspot/coldspot results
within 1e-10; default FTICR predictions remain equal. Deliberate boundary and
mean differences have separate tests. Both continuous and processed imzML,
memory modes, backed storage, profile rejection, corrupt binaries, metadata,
annotations, and run/coordinate validation are exercised.

Controlled imzML fixtures explicitly use 64-bit intensities. CardinalIO's
default writer uses 32-bit intensities, which rounds these inputs before
extraction and causes small MPM density differences. Production import retains
the precision present in the file; no scientific comparison tolerance was
relaxed for this difference.

`BiocManager::valid()` returned TRUE. The package builds and installs. R CMD
check completed with zero errors, 10 warnings, and 2 notes. Existing regression
scripts, examples, and the new automated tests passed. Warnings concern legacy
filenames, non-ASCII plotting code, undeclared/unused imports, S3 method
signatures, documentation, data compression, and vignette packaging. Notes
include legacy global-function declarations; suggested rmarkdown was unavailable.
PDF manuals and vignette rebuilding were disabled for this check. Full results
are recorded in `package-check.txt` and the change log.
Run automated tests with `Rscript --vanilla docs/Cardinal-migration/run-tests.R`.

Initial local-file benchmark (not a large-scale acceptance claim):

| Measurement | Result |
|---|---:|
| Acquired pixels | 24,504 |
| Binary spectral file | 681,849,892 bytes |
| XML file | 59,735,434 bytes |
| File-backed import | 1.73 s |
| Gaussian extraction, m/z 788.5447 | 10.16 s |
| Repeated extraction | 10.35 s |
| Detected pixels | 17,646 |
| Cardinal object size | 7,180,208 bytes |
| Whole-process peak working set | 1,494,827,008 bytes |

Peak memory includes R, package/data loading, XML parsing, import, and both
extraction runs; it is not an incremental extraction memory estimate. During
the repeated extraction, Rprofmem recorded no individual R allocation >= 1 MB.
This does not measure native allocations or establish constant memory use.
Disk-cache state was uncontrolled, only one local dataset was measured, and
legacy import/extraction performance was not benchmarked. No scalability
improvement is claimed from these numbers. CPU/RAM/storage inventory could not
be retrieved through WMI in the restricted environment.

The raw benchmark file triggers a Cardinal warning about an unrecognized
scan-direction ontology term; metadata preservation is bounded by Cardinal's
reader behavior. The legacy processed example has per-spectrum centroid tags
without a file-level tag and is read with unknown centroid status. It is
rejected for extraction until that status is established explicitly. The
benchmark uses the original file with an explicit file-level centroid tag.
Raw local datasets were neither changed nor committed.

## Intermediary cleanup: explicit spatial dependencies and readable source

The release-era dynamic importer has been retired. Spatial calls now identify
their defining package directly; no scan of installed packages or copied
spatstat functions is needed when moleculaR loads. The direct dependencies are:

| Package | Responsibility in moleculaR |
|---|---|
| `spatstat.geom` | Point patterns, windows, images, geometry, and plotting |
| `spatstat.explore` | Spatial density estimation and image smoothing |
| `spatstat.random` | CSR reference sampling and simulated point patterns |
| `spatstat.utils` | Input checks and polygon utilities |

Ordinary generic calls such as `summary(window)` still use standard S3
dispatch. Direct method calls use their explicit owner. The scientific spatial
calculations and the saved legacy numerical fixtures are preserved. Frozen
legacy source is replayed only in a private reference environment that binds
its three historical geometry names.

The intermediary suite has **302 passing assertions across 18 tests**, including
the original numerical comparisons and new spatial/plotting coverage. The
installed child-package versions and environment validity are recorded in
`environment.txt`; the retired import and meta-package entries are removed from
`renv.lock`.
The source archive also builds and installs successfully. R CMD check passes
all 302 assertions and the older standalone regression scripts with zero errors;
the 10 remaining warnings and 2 notes are recorded in `package-check.txt`.

For source navigation, start with these files:

| File | Purpose |
|---|---|
| `R/readMsi.R` | Public import: validate the file pair, read Cardinal data, select a run |
| `R/msi-accessors.R` | Read metadata, select pixels, validate the spatial/processing state |
| `R/msi-extractIon.R` | Validate targets, aggregate spectra, and construct point patterns |
| `R/searchAnalyte.R` | Public dispatch between Cardinal and the temporary legacy path |
| `tests/testthat/helper-cardinal.R` | Small synthetic inputs and frozen-reference adapters |
| `tests/testthat/test-spatstat-namespaces.R` | Coverage of spatial operations without the importer |

The helpers contain explanatory comments on coordinate order, mean denominators,
inclusive mass boundaries, annotation matching, processing queues, and memory
access. Development scripts begin with purpose, invocation, and output notes.
From the repository root, use `Rscript --vanilla` followed by the script path:

| Script under `docs/Cardinal-migration/` | When to use it |
|---|---|
| `setup-environment.R` | Install compatible development dependencies into the local library |
| `probe-environment.R` | Inspect installed versions and generic signatures without changing them |
| `inspect-api.R` | Inspect method bodies and exercise two tiny native spectra |
| `run-tests.R` | Run automatic regressions and update the per-test CSV |
| `document-package.R` | Regenerate API help, namespace, and source collation |
| `check-package.R` | Build/check a source archive and save diagnostic reports |
| `validate-environment.R` | Verify the Bioconductor environment and intentionally update renv.lock |
| `benchmark-first-milestone.R` | Measure an existing local file-backed dataset |
| `generate-legacy-fixtures.R` | Deliberately rebuild the historical reference; never a routine test step |

## Remaining migration

The single-analyte milestone and its original validation above are historical.
Subsequent work is described in [phase-3.md](phase-3.md) for batch/lipid
extraction and [phase-4.md](phase-4.md) for measured FWHM and preprocessing.
Shiny integration is described in [phase-5.md](phase-5.md). Example-data cleanup, removal of legacy dependencies,
cross-platform testing, and large-scale validation remain later phases.
