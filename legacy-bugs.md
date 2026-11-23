# Legacy bug register

Track defects discovered while converting moleculaR from MALDIquant to Cardinal.
This is the repair backlog for a **separate legacy working copy**. A fix on
`Cardinal-dev` does not mean that the corresponding bug has been fixed on `main`.

Keep stable IDs when adding findings. Record suspected bugs as unconfirmed until
there is a reproducible failure or sufficient evidence. Intentional changes in
preprocessing or data representation belong in `scientific-decisions.md`; they
are not automatically legacy bugs. Package-check housekeeping is recorded in
the change log and check reports.

## Verified legacy baseline

- Date verified: 2026-11-23.
- Local `main` commit: `a2d7f3500da5928275eb651f5c86b33beb189d67`.
- Migration reference commit: `ce10618b608e840f888bb71bb68794c8b44535ab`,
  preserved on `legacy/MALDIquant-reference`.
- Source was inspected with `git show main:<path>` and replayed in private
  environments. `main` was neither checked out nor modified. Remote branch
  freshness and protection settings were not changed or independently checked.
- Reproduction environment: R 4.6.1, MALDIquant 1.22.3, spatstat.geom 3.8.3.
  See `environment.txt` for the wider migration environment.

Verify the actual commit in the future legacy checkout before applying a fix.
To read this register from another checkout of the same repository:

```text
git show Cardinal-dev:docs/Cardinal-migration/legacy-bugs.md
```

## Backlog

| ID | Defect | Priority | Legacy status | Cardinal-dev status |
|---|---|---|---|---|
| LEGACY-001 | Nearest-index selection admits masses outside the intended search interval | High: scientific correctness | Fixed and validated on main, 2026-11-23; see repair log below | Fixed for Cardinal input in `321e1bc`; temporary sparse path still affected |
| LEGACY-002 | Simulation window operations depend on unimported functions being attached | Medium: affected simulation shapes fail | Fixed and validated on main, 2026-11-23; see repair log below | Fixed in `e111d98` |
| LEGACY-003 | A valid rectangular `spwin` causes an error in searchAnalyte | Medium: supported spatial input fails | Confirmed; open on verified main | Handled by Cardinal input in `321e1bc`; temporary sparse path still affected |

## LEGACY-001 — Mass selection can extend beyond three sigma

- **Discovered:** 2026-11-23, first-milestone legacy/Cardinal numerical comparison.
- **Affected code:** `R/searchAnalyte.R`; the `idx`, `idxlwr`, `idxupr`,
  `idxlwr:idxupr`, and corresponding Gaussian-weight indexing expressions.
- **Affected entry points:** single-analyte search and callers such as batch and
  lipid searches that use this function.
- **Legacy repair commit / validation:** `006dc71` (`fix: LEGACY-001 enforce inclusive three-sigma mass selection`); validated 2026-11-23 (details below).

### Cause and expected behavior

The intended interval is inclusive:

```r
sigma <- fwhm / (2 * sqrt(2 * log(2)))
tolerance <- 3 * sigma
selected <- which(abs(spData$mzAxis - m) <= tolerance)
```

Legacy code instead finds the global-axis features nearest the two interval
endpoints and takes the complete span between those feature indices. A nearest
endpoint feature can be outside the intended interval, introducing extra peaks.
The new Cardinal extraction path uses the numerical condition above directly.
This rule is moleculaR's scientific definition, not a restriction imposed by
Cardinal itself.

### Minimal reproduction in a legacy installation

Use a fresh R session with the legacy moleculaR package installed. The polygonal
window below deliberately avoids LEGACY-003 so that mass selection is isolated.

```r
library(moleculaR)

measuredMz <- c(99.90, 100.00, 100.14)
measuredIntensity <- c(2, 10, 8)

spData <- sparseIntensityMatrix(
    spmat = Matrix::Matrix(matrix(measuredIntensity, nrow = 1), sparse = TRUE),
    mzAxis = measuredMz,
    coordinates = data.frame(x = 1, y = 1)
)

rectangle <- spatstat.geom::owin(c(0.5, 1.5), c(0.5, 1.5))
polygon <- spatstat.geom::as.polygonal(rectangle)

result <- suppressWarnings(searchAnalyte(
    m = 100,
    fwhm = 0.1,
    spData = spData,
    wMethod = "sum",
    spwin = polygon
))

result$marks$intensity
# Verified legacy result: 20
# Correct strict-window result: 12
```

The interval is approximately `[99.8726, 100.1274]`. The peak at `100.14`
contributes intensity 8 in the legacy result despite lying outside that interval.

### Scientific impact

- `sum` includes outside intensities directly; `max` can select an outside peak.
- `mean` includes extra feature columns in both the sum and the denominator,
  so the result can increase or decrease.
- Gaussian weighting attenuates the error but does not eliminate it. The weight
  at exactly three sigma is `exp(-4.5)`, approximately **0.011109**, or **1.11% of
  that peak's intensity**. An outside peak has a lower weight, but a sufficiently
  intense peak can still contribute substantially. This is not a bound on the
  relative error in the final analyte intensity.
- Detected pixels and downstream maps can be affected. The size of that effect
  on historical biological datasets has not been quantified. This finding alone
  does not establish that every previous map is materially different.

### Legacy fix and acceptance criteria

1. Select global feature columns using the inclusive numerical interval.
   Use the same selected indices for intensities and Gaussian weights.
2. Preserve the Gaussian formula and other spatial/annotation behavior. Preserve
   the legacy shared-feature interpretation of `mean`; this bug fix should not
   introduce the unaligned-array mean definition as an unrelated change.
3. Preserve the existing no-match return convention when no column is selected.
4. Add regressions for exact endpoints, immediately outside both endpoints,
   coarse/nonuniform axes, no matches, and all four aggregation methods.
5. Include an unusually intense outside peak to verify that Gaussian attenuation
   cannot conceal accidental inclusion. Keep zero-valued acquired pixels in the
   tissue domain.
6. Reassess representative legacy maps with fixed seeds and document any changes
   in extracted intensities, detected pixels, densities, and classifications.

**Existing migration evidence:**
`tests/testthat/test-cardinal-legacy-spatial.R`, test
`strict mass windows quantify the nearest-index boundary difference`.
Retain the original legacy fixture as historical evidence; add corrected legacy
expectations separately when fixing this defect.

### Legacy repair log (2026-11-23)

Verified this checkout was on `main` at the recorded `a2d7f35` baseline before
editing. Reproduced the original result of **20** in a fresh installed-package
session with only moleculaR attached; the corrected result is **12**.

`searchAnalyte()` now selects feature columns using
`which(abs(spData$mzAxis - m) <= threeS)` and uses those same columns for Gaussian
weights. The Gaussian formula, shared-feature denominator for `mean`,
confirmation/metadata behavior, spatial windows, and empty-`ppp` return convention
are preserved. The search interval is also documented in the function help.

Validation used R 4.6.1, Matrix 1.7.5, MALDIquant 1.22.3, import 1.3.4, and
spatstat.geom 3.8-3. Dependencies were read from the existing local library;
the package was installed into an isolated validation library. The pre-existing
`renv.lock` changes were left untouched. No spatstat package was attached.

- `Rscript --vanilla tests/legacy001-mass-window.R`: **passed** after failing
  against the baseline. Covers exact endpoints, adjacent outside doubles on
  both sides, outside intensities of `1e12`, coarse/nonuniform axes, singleton
  selection, no matches/empty axes/all-zero signal, all four aggregation methods,
  zero-signal acquired pixels in the tissue window, metadata, confirmation, and
  batch search.
- Existing `tests/useDefaultFwhm.R`, `tests/estimateFwhm.R`,
  `tests/ionImage2spp.R`, and `tests/shinyFwhm.R`: **all passed**. Existing
  point-rejection warnings from probability-map masks and a Shiny package-build
  version warning were observed.
- `tests/manual/legacy001-map-comparison.R` compares the installed correction
  with the original function read from the immutable baseline via `git show`.
  It uses the bundled **24,504 acquired pixels**, the bundled measured FWHM
  model, unchanged tissue windows, **seed 1**, **bandwidth 2**, and all four
  aggregation methods. Full numeric output is retained in
  `tests/manual/legacy001-map-results.csv`.

Representative Gaussian-map results (before -> after):

| m/z | Detected pixels | Changed intensity pixels | Maximum intensity change | Maximum analyte-density change | Hotspot pixels | Coldspot pixels | Changed classifications |
|---|---|---|---|---|---|---|---|
| 788.5447 | 17643 -> 17643 | 0 | 0 | 0 | 5829 -> 5829 | 7343 -> 7343 | 0 |
| 760.5850422936694 | 24262 -> 24261 | 113 | 62665.1383 | 1.69605e-9 | 2603 -> 2621 | 2739 -> 2874 | 153 |
| 496.33972664966376 | 24441 -> 24434 | 7 | 4723.9971 | 2.24491e-9 | 4464 -> 4590 | 5672 -> 5952 | 406 |

All four methods leave the 788.5447 example unchanged. For 760.58504, the
`sum`/`mean`/`max` maps change 133/133/144 classifications; for 496.33973,
they each change 391. Removing outside columns changes the shared-feature mean
denominator, so many mean intensities legitimately increase. Changing detection
counts also changes CSR sampling even under the same seed; maximum CSR-density
changes for the two Gaussian maps are 3.55461e-5 and 6.08161e-5. Thus small
weighted-intensity changes need not imply identical classifications. These
results describe this bundled dataset and these mapping parameters.

## LEGACY-002 — Two simulation shapes use missing function lookups

- **Discovered:** 2026-11-23, intermediary spatstat namespace audit.
- **Affected code:** `R/simulation-functions.R`, `createHotspotWin()`;
  `R/manualSpatstatImport.R` omits the two required operations.
- **Legacy repair commit / validation:** `fix: LEGACY-002 import simulation window operations dynamically`; validated 2026-11-23 (details below). This entry is included in the repair commit.

### Reproduction and cause

In a fresh session, attach only the legacy moleculaR package:

```r
library(moleculaR)

createHotspotWin("multiCircle", radius = 2, backgroundLength = 20)
# Verified failure: could not find function "union.owin"

createHotspotWin("donut", radius = c(1, 2), backgroundLength = 20)
# Verified failure: could not find function "setminus.owin"
```

The function passes those names as strings to `do.call()`, but the manual
importer does not import either function. Attaching an appropriate spatstat
package can mask the problem, making behavior depend on the user's session.
The failures were reproduced from the verified main source with the geometry
binding supplied by its historical importer and without attaching spatstat.

### Legacy fix and acceptance criteria

- Add `union.owin` and `setminus.owin` to the existing dynamic spatstat importer
  and pass those imported function objects to `do.call()`. This follows the
  legacy branch's cross-version import policy; it supersedes the original
  migration-specific suggestion to hard-code a spatstat subpackage namespace.
- Test `centralCircle`, `multiCircle`, and `donut` in a fresh session that does
  not attach spatstat. Verify geometry and ensure simulation callers work.
- Keep shape definitions, radii, and spatial units unchanged.

**Existing migration evidence:** commit `e111d98` and
`tests/testthat/test-spatstat-namespaces.R`, test
`all simulation shapes resolve their geometry explicitly`.

### Legacy repair log (2026-11-23)

Reproduced both missing-function errors with the installed baseline and no
spatstat packages attached. Added `union.owin` and `setminus.owin` to
`spatstatFuncs` in `manualSpatstatImport.R`. `createHotspotWin()` passes those
imported functions directly to `do.call()`. The existing `import` and `spatstat`
dependencies suffice; no subpackage namespace or new dependency is introduced.

Reinstalled the package into the isolated validation library and ran
`Rscript --vanilla tests/legacy002-simulation-imports.R`: **passed**, after
reproducing a missing `union.owin` error before the correction. The test asserts
that no spatstat package is attached at the beginning or end. It checks:

- Central circle, five-circle union, and donut areas/bounds/units, including the
  donut hole, reversed donut radii, distinct circle radii, and overlapping circles.
  At radius 2 and background length 20, areas are approximately 12.56132463,
  62.80662308, and 9.42099345 (donut radii 1 and 2). Polygon-clipping comparisons
  allow relative tolerance `1e-7`.
- `simulatespp()` for all three shapes with seed 1: deterministic coordinates
  and intensities, foreground/background intensity quantiles, metadata, and
  unchanged background area of 400.
- The same simulation callers using measured intensities from the first 100
  spectra of the bundled MSI dataset; returned intensities come from that data.

The LEGACY-001 regression also passed after the importer change. Shape
definitions, radii, units, sampling, and intensity assignment are unchanged.

## LEGACY-003 — Rectangular spatial windows fail in legacy search

- **Discovered:** 2026-11-23, first-milestone boundary regression setup.
- **Affected code:** `R/searchAnalyte.R`, the supplied-window diagnostic beginning
  with `owinInfo <- summary(spwin)`.
- **Legacy repair commit / validation:** not yet supplied.

### Reproduction and cause

Reuse `spData` and `rectangle` from LEGACY-001:

```r
searchAnalyte(m = 100, fwhm = 0.1, spData = spData, spwin = rectangle)
# Verified legacy failure: argument is of length zero
```

A rectangular `owin` is valid under the documented input contract, but its
summary does not supply the polygon fields assumed by the diagnostic.
`if (owinInfo$npoly > 1)` consequently receives a length-zero condition.
Supplying the equivalent polygon avoids this error; it should not be required
merely to use an otherwise valid spatial window.

### Legacy fix and acceptance criteria

- Validate the `owin` class and either normalize the window to polygonal form
  before polygon-specific diagnostics or restrict those diagnostics to polygonal
  windows. Preserve the spatial domain.
- Test a rectangle, an equivalent polygon, a disconnected/holey polygon, invalid
  input, and a valid no-match search. Rectangular and equivalent polygonal
  windows should give the same meaningful scientific result.
- Keep this repair separate from the mass-selection correction so each can be
  tested independently.

**Existing migration evidence:** Cardinal input normalizes a valid window with
`spatstat.geom::as.polygonal()` in `R/msi-extractIon.R`. Migration annotation and
no-match tests already use rectangular windows successfully.

## Maintaining the register

For each new entry record discovery date/context, affected functions and verified
legacy revision, confirmation status, minimal reproduction, observed/expected
behavior, cause, scientific or operational impact, suggested repair, acceptance
tests, and independent legacy/Cardinal repair commits.

When repairing the legacy package in another working copy, update the **legacy**
status only after the corrected tests run there. Retain the original reproduction
and record numerical changes. Link the eventual legacy fix commit and, when
applicable, the commit merged into `main`. Do not close an entry solely because
the Cardinal implementation passes its tests.
