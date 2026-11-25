# Dependency and data representation audit

Starting commit: `ce10618b608e840f888bb71bb68794c8b44535ab` on `Cardinal-dev`.
Preserved as `legacy/MALDIquant-reference`. Audit covers R, inst, vignettes,
tests, DESCRIPTION, NAMESPACE, and repository-wide searches for peak processing.

| Area | Dependencies and assumptions | Migration action |
|---|---|---|
| DESCRIPTION / NAMESPACE | MALDIquant imported; MALDIquantForeign suggested; explicit Collate; generated exports | Add Cardinal and APIs directly called; preserve old dependencies until validated cleanup |
| readCentrData | XML first-50-line continuous check; importImzMl returns MassPeaks | Add readMsi; retain old reader for now |
| readSingleSpect | Two-column TSV becomes MassSpectrum | Numeric input rewrite in phase 4 |
| createSparseMat / S3classes | All masses/intensities unlisted; global union axis; pixels by features; Matrix dgCMatrix | Bypass for Cardinal; retain reference utility |
| searchAnalyte | match.closest for confirmation and inclusive nearest-index span; dense selected columns; global mean denominator | Add strict Cardinal extraction; document boundary and mean differences |
| batchAnalyteSearch | Temporary sparse path still uses coordinates/mzAxis and mclapply | Phase 3 adds native Cardinal pixel chunks, sparse query results, BiocParallel and optional cache |
| batchLipidSearch / legacy-searchLipid | Sparse assumptions remain on the temporary path and archived searchLipid | Phase 3 routes native lipid/adduct queries through the shared batch engine and preserves metadata |
| normalizeIntensity | MALDIquant intensity setters; mergeMassPeaks; exact shared feature matching; IS imputation | Phase 4 validate factors/scaling and missing-standard policy |
| filterPeaks | Creates sparse matrix; stored-entry frequencies via diff(@p); keeps >= threshold; reconstructs MassPeaks/FWHM metadata | Aligned Cardinal filtering in phase 4; code threshold differs from prose <= statement |
| estimateFwhm | MassSpectrum/MassPeaks tests/accessors, SuperSmoother peak detection, match.closest, numeric half-height interpolation | Preserve fitting/defaults; replace adapters and validate peak selection in phase 4 |
| useDefaultFwhm / getFwhm | FTICR interpolation loads bundled processed-example-Data; no peak processing in default path | Preserve numerical reference independently before data cleanup |
| gaussWeight | Base Gaussian density ratio | Preserve unchanged |
| createSpatialWindow | x/y pixel grid -> polygonal owin; optional explicit pixelSize; no spectral dependency | Reuse with all acquired coordinates |
| ionImage2spp | Numeric image/array and optional msiSlices centre attribute; no runtime MALDIquant calls | Preserve; MALDIquant examples/tests retained for now |
| probMap / superImposeAnalytes / subsetAnalytes | ppp marks and metadata identifier linkage, spatial window; no spectral objects | Preserve unchanged and compare scientific fields with fixed seeds |
| plotting-functions / superimposeNoise | MALDIquant plotting; sparse axes; reconstruction of MassPeaks | Later migration needed, outside first milestone |
| recreateExData | Legacy import/binning/filtering/sparse conversion and .shrinkData accessors | Phase 6; preserve FTICR data first |
| runShiny / inst/shiny | Launcher independent; apps build sparse data, snap targets to mzAxis, use MALDIquant coordinates | Phase 5; package-app already copies paired uploads into session storage, preserve that work |
| vignettes | Walkthrough Rmd, R and .orig demonstrate legacy objects and pipeline | Phase 6 update active workflows |
| tests | Interactive/manual FWHM scripts and standalone ionImage2spp regression; no testthat harness | Add dedicated automatic Cardinal tests and immutable legacy fixtures |
| classDefinitions | Commented proposed S4 classes only | No runtime dependency |

Cardinal shared-axis matrices are **features by pixels**, the transpose of the
legacy sparse representation. Arrays have individual mass axes and must not be
forced through global alignment. File-backed arrays must remain authoritative;
never unlist the complete dataset or convert its full intensity matrix.

Additional audit hits in man/README and bundled data are documentation or legacy
fixtures; they still require the phase 6 cleanup audit. Local `.package.references`
contains supplied documentation, not package runtime code. Raw external imzML/ibd
pairs are available in ignored `data-external/`; synthetic fixtures and generated
small imzML pairs provide redistributable first-milestone coverage.

## Intermediary update after the first milestone

Removed `R/manualSpatstatImport.R`, its documentation topic, and its Collate /
roxygen include references. Runtime `import::here()` and the `spatstat`
meta-package dependency are replaced by explicit child-package calls and direct
Imports for spatstat.geom, spatstat.explore, spatstat.random, and spatstat.utils.
The qualification pass updated 94 parser-identified references without changing
their normalized parsed expressions. Two additional string-based simulation
lookups now pass `spatstat.geom::union.owin` and `spatstat.geom::setminus.owin`
directly to `do.call()`. Existing explicitly namespaced calls are preserved.

The legacy source fixture remains unchanged. Its private replay environment
provides only ppp, owin, and as.polygonal; production namespaces no longer contain
copied spatstat functions or installed-package enumeration objects.

## Phase 3 update (2026-11-25)

`R/msi-extractIons.R` owns multi-query extraction, query validation and
deduplication, sparse result assembly, optional caching, pixel chunk sizing,
BiocParallel backend ownership, common spatial preparation and lipid query
metadata. Both public batch functions dispatch native MSI objects there; legacy
sparse implementations remain intact. No Cardinal dataset is converted into a
global spectral matrix. The result matrix contains extracted query values only.

The cache reads backing-file paths through `matter::path()`, so matter moves
from Suggests to Imports with the already-validated version floor. No package
version or lockfile changes are required. Spectral metadata use the exported
`ProtGenerics::spectraData()` accessor and integer indexing; coercing the entire
S4 metadata table to a list would invoke an unnecessary optional IRanges path.
New source is included in Collate and batch help is generated through roxygen2.
See `phase-3.md`, the new batch regression file, and `batch-benchmark.csv` for
the validation and measured scope.
