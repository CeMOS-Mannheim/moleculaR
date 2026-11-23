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
| LEGACY-001 | Nearest-index selection admits masses outside the intended search interval | High: scientific correctness | Confirmed; open on verified main | Fixed for Cardinal input in `321e1bc`; temporary sparse path still affected |
| LEGACY-002 | Simulation window operations depend on unimported functions being attached | Medium: affected simulation shapes fail | Confirmed; open on verified main | Fixed in `e111d98` |
| LEGACY-003 | A valid rectangular `spwin` causes an error in searchAnalyte | Medium: supported spatial input fails | Confirmed; open on verified main | Handled by Cardinal input in `321e1bc`; temporary sparse path still affected |

## LEGACY-001 — Mass selection can extend beyond three sigma

- **Discovered:** 2026-11-23, first-milestone legacy/Cardinal numerical comparison.
- **Affected code:** `R/searchAnalyte.R`; the `idx`, `idxlwr`, `idxupr`,
  `idxlwr:idxupr`, and corresponding Gaussian-weight indexing expressions.
- **Affected entry points:** single-analyte search and callers such as batch and
  lipid searches that use this function.
- **Legacy repair commit / validation:** not yet supplied.

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

## LEGACY-002 — Two simulation shapes use missing function lookups

- **Discovered:** 2026-11-23, intermediary spatstat namespace audit.
- **Affected code:** `R/simulation-functions.R`, `createHotspotWin()`;
  `R/manualSpatstatImport.R` omits the two required operations.
- **Legacy repair commit / validation:** not yet supplied.

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

- Pass `spatstat.geom::union.owin` and `spatstat.geom::setminus.owin` as function
  objects to `do.call()`, with the appropriate direct dependency declared.
- Test `centralCircle`, `multiCircle`, and `donut` in a fresh session that does
  not attach spatstat. Verify geometry and ensure simulation callers work.
- Keep shape definitions, radii, and spatial units unchanged.

**Existing migration evidence:** commit `e111d98` and
`tests/testthat/test-spatstat-namespaces.R`, test
`all simulation shapes resolve their geometry explicitly`.

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
