# moleculaR: Migration from MALDIquant to Cardinal

## Comprehensive Development Specification and Instructions for an AI Coding Agent

**Project:** moleculaR  
**Repository:** CeMOS-Mannheim/moleculaR  
**Development branch:** `Cardinal-dev`  
**Migration strategy:** Complete replacement of the MALDIquant-based MSI infrastructure with Cardinal  
**Legacy strategy:** Preserve the MALDIquant-based implementation in a separate Git branch  
**Primary language:** R  
**Target infrastructure:** Cardinal 3, matter, Bioconductor  
**Document purpose:** Autonomous development specification, implementation roadmap, and validation protocol

---

# 1. Executive Summary

## 1.1 Project objective

The objective of this development project is to migrate the moleculaR R package from its existing MALDIquant-based mass spectrometry imaging (MSI) infrastructure to a modern Cardinal-based architecture.

The migration must preserve the scientific functionality of moleculaR while substantially improving its ability to process large MSI datasets.

The current package relies on MALDIquant and MALDIquantForeign for:

- Importing MSI datasets.
- Representing individual mass spectra.
- Spectral preprocessing.
- Peak detection and alignment.
- Peak intensity normalization.
- Estimating spectral peak widths.
- Extracting analyte intensities.

The package subsequently converts the MALDIquant representation into a custom sparse intensity matrix using the Matrix package.

This architecture works well for relatively small or moderate MSI datasets but becomes inefficient when dealing with large datasets.

Cardinal provides a more suitable infrastructure through:

- Native imzML reading and writing.
- File-backed spectral storage using matter.
- Efficient access to large spectral datasets.
- S4 classes specifically designed for MSI experiments.
- Deferred preprocessing.
- Chunk-based computation.
- Peak detection and alignment.
- Spectral normalization.
- Parallel processing through BiocParallel.
- Preservation of MSI-specific metadata.

The migration must use these capabilities directly.

**The objective is not simply to replace MALDIquant function calls with corresponding Cardinal functions.**

Instead, the objective is to redesign moleculaR's spectral data-handling architecture around Cardinal's native data structures and computational model.

## 1.2 Selected migration strategy

The selected strategy is a complete migration to Cardinal.

MALDIquant will not remain a permanently supported alternative backend in the main development branch.

However, the existing implementation must be preserved and used as a scientific reference during migration.

During development, MALDIquant may be installed in a separate legacy environment for numerical comparisons.

The final Cardinal-based package should not require MALDIquant or MALDIquantForeign for its normal operation.

## 1.3 Scientific preservation requirement

The molecular probability mapping methodology implemented by moleculaR must remain scientifically unchanged unless a modification is explicitly required to correct an identified issue.

The following functionality must be preserved:

- Single-analyte detection.
- FWHM-dependent mass matching.
- Gaussian-weighted analyte intensity extraction.
- Alternative intensity aggregation methods.
- Molecular probability map generation.
- Kernel density estimation.
- Complete spatial randomness reference models.
- Hotspot and coldspot detection.
- Spatial statistical testing.
- False discovery rate correction.
- Collective projection probability maps.
- Batch analyte searches.
- SwissLipids-based searches.
- METASPACE-based annotation verification.
- Spatial visualization.
- Shiny-based workflows.

Changes to spectral preprocessing may produce differences between the old and new implementations.

Such differences must be identified, quantified, and documented.

They must never be introduced silently.

---

# 2. Instructions for the Development Agent

## 2.1 General operating principles

The development agent must follow these principles throughout the migration.

1. Work exclusively on the `Cardinal-dev` branch.
2. Inspect the existing source code before modifying any function.
3. Preserve the legacy implementation as a reproducible reference.
4. Make incremental, logically separated changes.
5. Add tests for every modified or newly implemented function.
6. Preserve existing public APIs where feasible.
7. Avoid unnecessary changes to the scientific algorithms.
8. Use Cardinal's native data structures.
9. Avoid materializing entire MSI datasets unless explicitly justified.
10. Prefer file-backed and chunk-based operations.
11. Maintain compatibility with Linux, Windows, and macOS where supported by Cardinal.
12. Record all scientific and technical decisions.
13. Do not remove MALDIquant dependencies until the corresponding functionality has been replaced and validated.
14. Do not declare the migration complete until the acceptance criteria in this document are satisfied.
15. When suitable alternatives exist, prioritize Bioconductor packages over CRAN packages. Verify compatibility and scientific suitability before selecting a dependency.
16. Document every development step in `docs/Cardinal-migration/change-log.md`, including files and functions changed, validation results, numerical differences, limitations, and remaining work. Update this log as each step is completed.

## 2.2 Required development workflow

For each major implementation stage:

1. Inspect the relevant existing moleculaR functions.
2. Identify their current behavior and dependencies.
3. Examine the corresponding Cardinal functionality.
4. Define the desired replacement behavior.
5. Implement the replacement.
6. Write unit tests.
7. Compare results against the legacy implementation.
8. Evaluate memory consumption and computational performance.
9. Update documentation.
10. Commit the completed stage.

Do not perform a wholesale rewrite of the entire package before testing the individual components.

## 2.3 Required reporting

After each development stage, provide a concise report containing:

- Files modified.
- Functions added or changed.
- Existing functionality preserved.
- Tests executed.
- Test results.
- Numerical differences from the legacy implementation.
- Performance observations.
- Known limitations.
- Remaining work.

Whenever an implementation decision could affect scientific results, document the decision before proceeding.

Use an explicit decision record rather than relying on comments scattered throughout the source code.

---

# 3. Source Code and Documentation to Inspect

## 3.1 Primary source repository

The primary source of truth is the actual moleculaR source code in the repository.

Inspect at least the following files:

| File | Main responsibility |
|---|---|
| `DESCRIPTION` | Package dependencies and compatibility |
| `NAMESPACE` | Exported functions and imports |
| `R/readCentrData.R` | imzML import |
| `R/readSingleSpect.R` | Import of individual profile spectra |
| `R/createSparseMat.R` | Construction of sparse MSI representation |
| `R/S3classes.R` | Custom object constructors |
| `R/searchAnalyte.R` | Single-analyte intensity extraction |
| `R/batchAnalyteSearch.R` | Batch analyte extraction |
| `R/batchLipidSearch.R` | SwissLipids searches |
| `R/estimateFwhm.R` | FWHM estimation |
| `R/useDefaultFwhm.R` | Instrument-dependent FWHM models |
| `R/gaussWeight.R` | Gaussian intensity weighting |
| `R/normalizeIntensity.R` | Normalization |
| `R/filterPeaks.R` | Peak-frequency filtering |
| `R/createSpatialWindow.R` | Spatial tissue boundaries |
| `R/ionImage2spp.R` | Ion-image to point-pattern conversion |
| `R/probMap.R` | Molecular probability mapping |
| `R/superImposeAnalytes.R` | Combining analytes |
| `R/subsetAnalytes.R` | Analyte subsetting |
| `R/recreateExData.R` | Example-data generation |
| `R/runShiny.R` | Shiny application entry point |
| `inst/shiny/package-app.R` | Main Shiny application |
| `inst/shiny/example-app.R` | Example Shiny application |

Inspect all other files for direct or indirect dependencies on MALDIquant objects.

Do not assume that a file is independent of MALDIquant merely because it does not explicitly import MALDIquant.

For example, a function may depend on the structure of an object originally constructed using MALDIquant.

## 3.2 Required dependency audit

Search the repository for at least:

```bash
grep -RIn "MALDIquant" R inst vignettes tests
grep -RIn "MassPeaks" R inst vignettes tests
grep -RIn "MassSpectrum" R inst vignettes tests
grep -RIn "sparseIntensityMatrix" R inst vignettes tests
grep -RIn "createSparseMat" R inst vignettes tests
grep -RIn "mzAxis" R inst vignettes tests
grep -RIn "spmat" R inst vignettes tests
grep -RIn "coordinates(" R inst vignettes tests
```

Also inspect:

```bash
grep -RIn "binPeaks" .
grep -RIn "detectPeaks" .
grep -RIn "match.closest" .
grep -RIn "mergeMassPeaks" .
grep -RIn "importImzMl" .
```

Record the findings in:

`docs/Cardinal-migration/dependency-audit.md`

The audit must identify both direct dependencies and assumptions about data representation.

## 3.3 Cardinal documentation

Use the official Cardinal documentation as the authoritative reference for Cardinal behavior.

Required materials:

- Cardinal reference manual.
- Cardinal 3 user guide.
- Cardinal documentation for MSI data structures.
- Cardinal documentation for preprocessing.
- matter documentation.
- BiocParallel documentation.
- CardinalIO documentation.

The supplied reference materials include:

- `Cardinal.pdf`, version 3.10.0.
- Segmentation of MS imaging experiments with Cardinal.
- Classification of MS imaging experiments with Cardinal.
- Class comparison for MSI experiments with Cardinal.
- Cardinal and matter source material.

These references are useful for understanding the architecture and supported APIs.

**Important:** The supplied Cardinal manual documents version 3.10.0. At the time this development specification was prepared, the Bioconductor release package page identified Cardinal 3.14.0. The agent must verify the installed versions and the actual method signatures before implementation.

Do not mix assumptions from different Cardinal versions.

Record the selected versions in the package documentation and dependency lockfile.

---

# 4. Current moleculaR Architecture

## 4.1 Existing processing pipeline

The present implementation follows this general workflow:

```text
imzML + ibd
     |
     v
readCentrData()
     |
     v
MALDIquant::MassPeaks
     |
     v
MALDIquant::binPeaks()
     |
     v
Optional peak filtering
     |
     v
Optional normalization
     |
     v
createSparseMat()
     |
     v
sparseIntensityMatrix
     |
     v
searchAnalyte()
     |
     v
analytePointPattern
     |
     v
probMap()
     |
     v
molProbMap
     |
     v
Visualization / spatial analysis
```

The key problem is that the processing pipeline requires a complete MALDIquant representation before constructing the global sparse matrix.

For very large datasets, both operations can become computationally expensive.

## 4.2 Existing sparse matrix representation

The custom `sparseIntensityMatrix` class contains:

```r
list(
    spmat = ...,
    mzAxis = ...,
    coordinates = ...
)
```

Its intensity matrix has the orientation:

```text
              m/z 1    m/z 2    m/z 3    ...    m/z N

Pixel 1         I        I        I                I
Pixel 2         I        I        I                I
Pixel 3         I        I        I                I
...
Pixel M         I        I        I                I
```

The matrix is created by collecting the peak masses and intensities from all spectra.

The implementation relies on operations equivalent to:

```r
mass <- unlist(lapply(x, MALDIquant::mass))

intensity <- unlist(lapply(x, MALDIquant::intensity))

uniqueMass <- sort(unique(mass))
```

and subsequently constructs a `Matrix::sparseMatrix`.

This creates several limitations:

- Large intermediate vectors are generated.
- The entire dataset must be traversed.
- Many individual R objects are retained during preprocessing.
- A global m/z axis must be constructed.
- Global peak alignment is effectively required.
- The initial conversion may consume substantial memory.

Although the resulting matrix is sparse, creating it is not necessarily memory-efficient.

## 4.3 Existing analyte extraction

The current `searchAnalyte()` requires a `sparseIntensityMatrix`.

It:

1. Receives a theoretical analyte mass.
2. Receives an FWHM value at that mass.
3. Calculates a Gaussian standard deviation.
4. Determines the relevant mass range.
5. Identifies corresponding columns in the sparse matrix.
6. Extracts the spectral intensities.
7. Applies the requested intensity aggregation method.
8. Identifies pixels with positive analyte intensity.
9. Converts those pixels into an `analytePointPattern`.
10. Adds verification and analyte metadata.

This function is the most important target for migration.

## 4.4 Existing spatial analysis

The `probMap()` function accepts a spatial point pattern.

It does not fundamentally require a MALDIquant spectrum list.

This creates a useful separation between:

- Spectral data handling.
- Molecular probability mapping.

**The migration must preserve this separation.**

The spatial probability mapping code should continue to operate on existing moleculaR spatial objects.

The MSI data backend should not be introduced into the internal calculations of `probMap()` unless strictly necessary.

---

# 5. Target Architecture

The desired architecture is:

```text
                  imzML / ibd
                       |
                       v
             Cardinal::readMSIData()
                       |
                       v
                 Cardinal object
                  file-backed
                       |
            +----------+----------+
            |                     |
            v                     v
   MSImagingArrays       MSImagingExperiment
   individual m/z axes   shared m/z axis
            |                     |
            +----------+----------+
                       |
                       v
              moleculaR MSI access
                       |
            +----------+----------+
            |                     |
            v                     v
    Single-analyte           Batch-analyte
       extraction              extraction
            |                     |
            +----------+----------+
                       |
                       v
                Ion intensities
                       |
                       v
              analytePointPattern
                       |
                       v
                    probMap()
                       |
                       v
                   molProbMap
                       |
                       v
           Visualization / CPPMs
```

## 5.1 Architectural principles

The main Cardinal object must remain the authoritative representation of the MSI experiment.

Do not immediately convert it into a custom moleculaR matrix.

Avoid duplicating:

- m/z information.
- Pixel coordinates.
- Run information.
- Spectrum metadata.
- Spectral intensity arrays.

Use Cardinal accessors to retrieve the relevant information.

The custom moleculaR representation should be limited to information needed for the package's own analyses.

## 5.2 Preferred data representation

Use `MSImagingArrays` as the initial preferred representation for unaligned centroided data.

Support `MSImagingExperiment` for aligned datasets.

Both classes should be accepted by the primary analyte-extraction interface.

Do not require users to convert `MSImagingArrays` into `MSImagingExperiment` before generating molecular probability maps.

## 5.3 Proposed internal interface

Implement a small set of internal helpers.

Suggested function names:

```r
.msiCoordinates()
.msiRun()
.msiMassRange()
.msiNumberOfPixels()
.msiIsCentroided()
.msiExtractIon()
.msiExtractIons()
.msiSelectPixels()
.msiValidate()
```

These are proposed moleculaR functions, not existing Cardinal functions.

The public API should remain compact.

Avoid introducing unnecessary wrapper classes.

For Cardinal objects, the internal helpers should use Cardinal's S4 accessors.

---

# 6. Understanding Cardinal's Data Infrastructure

This section is essential. The development agent must understand these concepts before implementing the migration.

## 6.1 Cardinal and matter

Cardinal provides the high-level MSI analysis interface.

The matter package provides much of the underlying out-of-core data infrastructure.

The distinction is:

**Cardinal manages MSI-specific data structures and processing.**

**matter provides mechanisms for accessing and processing large arrays and vectors without loading everything into RAM.**

Cardinal combines these capabilities with MSI metadata and statistical functionality.

## 6.2 File-backed data

When an imzML dataset is read with:

```r
msData <- Cardinal::readMSIData(
    file = "dataset.imzML",
    memory = FALSE
)
```

the spectral contents can remain associated with the binary data on disk.

Cardinal creates data objects that reference the stored spectral data.

Only the required data are retrieved when an operation accesses them.

This differs substantially from reading an entire dataset into a list of ordinary R objects.

### Important limitation

File-backed storage does not mean every operation uses constant memory.

A function can still load a large amount of data if it requests a large subset or converts the complete object to an ordinary matrix.

The agent must inspect data-access patterns carefully.

Avoid operations such as:

```r
as.matrix(Cardinal::intensity(msData))
```

on large datasets.

Likewise, avoid converting the complete spectral data into a list of in-memory numeric vectors.

Such operations should only be used on explicitly small datasets or subsets.

## 6.3 The `MSImagingArrays` class

`MSImagingArrays` is designed for mass spectra where each spectrum may have its own m/z values.

This is particularly appropriate for centroided MSI datasets.

Conceptually:

```text
Pixel 1
    mz:        [100.1, 150.2, 200.4]
    intensity: [10.0, 20.0, 15.0]

Pixel 2
    mz:        [100.12, 175.3, 200.39]
    intensity: [12.0,  5.0, 18.0]

Pixel 3
    mz:        [100.11, 150.21, 300.5]
    intensity: [ 8.0,  9.0, 22.0]
```

The individual spectra do not need to contain identical m/z values.

This resembles MALDIquant's spectrum-list representation conceptually.

The important difference is that Cardinal can use file-backed storage.

### Construction

A small in-memory example:

```r
library(Cardinal)

mzList <- list(
    c(100.1, 150.2, 200.4),
    c(100.12, 175.3, 200.39),
    c(100.11, 150.21, 300.5)
)

intensityList <- list(
    c(10, 20, 15),
    c(12, 5, 18),
    c(8, 9, 22)
)

coordinates <- data.frame(
    x = c(1, 2, 3),
    y = c(1, 1, 1)
)

msData <- Cardinal::MSImagingArrays(
    spectraData = list(
        intensity = intensityList,
        mz = mzList
    ),
    pixelData = Cardinal::PositionDataFrame(coordinates),
    centroided = TRUE,
    continuous = FALSE
)
```

This type of object should be used for synthetic tests.

### Accessing spectral information

The main accessors are:

```r
Cardinal::spectraData(msData)
```

```r
Cardinal::mz(msData)
```

```r
Cardinal::intensity(msData)
```

```r
Cardinal::spectra(msData, "mz")
```

```r
Cardinal::spectra(msData, "intensity")
```

The returned arrays may be matter-backed.

A single spectrum can be accessed using the corresponding index.

For example:

```r
mzValues <- Cardinal::mz(msData)[[1]]

intensities <- Cardinal::intensity(msData)[[1]]
```

The index must refer to the same spectrum in both arrays.

### Pixel metadata

Use:

```r
Cardinal::pixelData(msData)
```

For coordinates:

```r
Cardinal::coord(msData)
```

For experimental runs:

```r
Cardinal::run(msData)
```

```r
Cardinal::runNames(msData)
```

### Important implications

The coordinate information and spectral information must remain synchronized.

When processing or subsetting the data, the agent must ensure that spectrum index `i` still corresponds to coordinate row `i`.

Never sort or reorder the coordinate table independently of the spectral data.

## 6.4 The `MSImagingExperiment` class

`MSImagingExperiment` is designed for datasets with a common m/z axis.

Conceptually:

```text
                  Pixel 1    Pixel 2    Pixel 3

m/z 100.1            10         12          8
m/z 150.2            20          0          9
m/z 200.4            15         18          0
m/z 300.5             0          0         22
```

The matrix orientation is:

**Rows = m/z features.**

**Columns = pixels/spectra.**

This is the transpose of the existing moleculaR sparse intensity matrix representation.

### Important consequence

Do not reuse indexing logic from the existing sparse matrix implementation without changing its orientation.

In the current moleculaR representation:

```r
spData$spmat[pixels, features]
```

In Cardinal's shared-axis representation:

```r
Cardinal::intensity(msData)[features, pixels]
```

### Accessing data

```r
mzAxis <- Cardinal::mz(msData)
```

```r
intensities <- Cardinal::intensity(msData)
```

```r
featureInfo <- Cardinal::featureData(msData)
```

```r
pixelInfo <- Cardinal::pixelData(msData)
```

The intensity array may be a matter-backed matrix.

### Subsetting

Features can be selected with:

```r
Cardinal::subsetFeatures(
    msData,
    mz >= 700,
    mz <= 900
)
```

Pixels can be selected using:

```r
Cardinal::subsetPixels(
    msData,
    x >= 10,
    x <= 100
)
```

These operations should generally be preferred over extracting complete matrices.

### Advantages

`MSImagingExperiment` is particularly useful for:

- Repeated feature extraction.
- Large predefined feature sets.
- Peak-frequency filtering.
- Feature-level statistics.
- Cardinal's segmentation and classification methods.
- Processing datasets whose peaks are already aligned.

### Limitations

Creating a common m/z axis requires alignment or resampling.

This may change the numerical representation of the original spectra.

For moleculaR, global alignment should not be mandatory.

## 6.5 Cardinal metadata

Cardinal distinguishes several types of metadata.

### `spectraData()`

Contains spectral arrays.

Examples:

- Intensities.
- m/z values.
- Additional spectral measurements.

### `featureData()`

Contains feature-level metadata.

Relevant primarily for `MSImagingExperiment`.

Examples:

- m/z values.
- Peak counts.
- Peak frequencies.
- Feature annotations.

### `pixelData()`

Contains spectrum- or pixel-level metadata.

Examples:

- x coordinate.
- y coordinate.
- Experimental run.
- Tissue annotation.
- Sample classification.
- Acquisition metadata.

### `experimentData()`

Contains experimental metadata associated with the MSI acquisition.

### `metadata()`

Contains additional general metadata.

### `processingData()`

Contains information about queued processing operations.

The agent must preserve these distinctions.

Do not discard experiment metadata when constructing processed datasets.

---

# 7. Cardinal Data Import

## 7.1 `readMSIData()`

The central import function is:

```r
Cardinal::readMSIData()
```

Recommended initial usage:

```r
msData <- Cardinal::readMSIData(
    file = "dataset.imzML",
    memory = FALSE,
    as = "MSImagingArrays"
)
```

This explicitly requests a representation suitable for spectra with individual m/z axes.

Using `as = "auto"` allows Cardinal to select the representation.

For moleculaR, explicit control is preferable during early development.

### Important arguments

| Argument | Purpose |
|---|---|
| `file` | Input MSI file |
| `memory` | Controls in-memory or file-backed access |
| `as` | Requested Cardinal representation |
| `check` | Optional binary-file verification |
| `mass.range` | Mass range used during conversion to a shared-axis representation |
| `resolution` | Shared-axis bin resolution |
| `units` | Absolute m/z or ppm units |
| `guess.max` | Number of spectra used when estimating a conversion axis |
| `parse.only` | Retrieve metadata without constructing the spectral object |
| `chunkopts` | Chunk-processing options |
| `BPPARAM` | Parallel processing configuration |

Some of these arguments are delegated to format-specific readers.

Verify the installed version's method signatures.

## 7.2 Continuous versus processed imzML

These concepts must not be confused with profile versus centroided spectra.

### Continuous imzML

Continuous imzML generally contains spectra sharing a common m/z axis.

Cardinal can represent such data as `MSImagingExperiment`.

### Processed imzML

Processed imzML permits different m/z arrays for different spectra.

Cardinal commonly represents such data as `MSImagingArrays`.

### Centroided data

Centroided data contain detected peaks rather than complete profile peak shapes.

### Profile data

Profile data contain sampled peak shapes.

**A processed imzML file is not necessarily centroided.**

This distinction is explicitly demonstrated in Cardinal's documentation.

### Required behavior

The initial Cardinal-native moleculaR analysis should support centroided MSI spectra.

For profile data:

- Do not silently treat the data as centroided.
- Validate the spectrum representation.
- Require or offer an explicit preprocessing workflow.
- Allow peak picking using Cardinal.

If centroid status is unknown, do not automatically assume the dataset is centroided.

Provide a documented validation mechanism.

---

# 8. Cardinal Spectral Processing

## 8.1 Deferred processing

Cardinal allows preprocessing steps to be queued.

For example:

```r
msData <- Cardinal::normalize(
    msData,
    method = "tic"
)
```

The operation is not necessarily executed immediately.

To execute queued processing:

```r
msData <- Cardinal::process(msData)
```

Multiple operations can be queued:

```r
msData <- msData |>
    Cardinal::normalize(method = "tic") |>
    Cardinal::smooth(method = "gaussian") |>
    Cardinal::reduceBaseline(method = "locmin") |>
    Cardinal::peakPick(SNR = 3)
```

These can subsequently be executed.

### Critical implementation requirement

moleculaR must define how it handles Cardinal objects containing queued processing operations.

Possible approaches:

1. Require users to process the data explicitly.
2. Process queued operations internally before analyte extraction.
3. Support deferred evaluation where the Cardinal API guarantees it.

**Recommended initial approach:**

Use explicit processing before molecular intensity extraction whenever queued preprocessing is present.

Do not silently call `reset()` or discard queued operations.

Record the processing history.

Be careful not to apply the same normalization or peak-picking operation twice.

## 8.2 Normalization

Cardinal supports:

```r
Cardinal::normalize(
    msData,
    method = "tic"
)
```

```r
Cardinal::normalize(
    msData,
    method = "rms"
)
```

```r
Cardinal::normalize(
    msData,
    method = "reference",
    ref = referenceMz
)
```

These correspond to:

- Total ion current normalization.
- Root-mean-square normalization.
- Reference-feature normalization.

They are not guaranteed to reproduce moleculaR's existing numerical behavior without testing.

Normalization factors and scaling conventions must be verified.

## 8.3 Peak picking

Cardinal provides:

```r
Cardinal::peakPick()
```

Example:

```r
msData <- Cardinal::peakPick(
    msData,
    method = "diff",
    SNR = 3
)

msData <- Cardinal::process(msData)
```

Supported peak-picking methods include:

- `diff`
- `sd`
- `mad`
- `quantile`
- `filter`
- `cwt`

These methods differ in their treatment of peak detection and noise estimation.

### Important instruction

Do not automatically perform peak picking on data that have already been centroided.

Doing so could alter the number, positions, or intensities of the detected peaks.

The initial Cardinal migration should preserve centroided input data unless the user explicitly requests additional preprocessing.

## 8.4 Peak alignment

Cardinal provides:

```r
Cardinal::peakAlign()
```

For example:

```r
aligned <- Cardinal::peakAlign(
    msData,
    tolerance = 200,
    units = "ppm"
)
```

This produces a shared-axis representation.

The resulting object is generally an `MSImagingExperiment`.

Peak alignment is useful when a common set of features is required.

However, its behavior differs from MALDIquant's `binPeaks()`.

The two functions must not be treated as numerically identical.

### Recommended approach

Make peak alignment optional.

Do not make it part of the default Cardinal-native `searchAnalyte()` workflow.

## 8.5 `peakProcess()`

Cardinal provides a combined workflow through:

```r
Cardinal::peakProcess()
```

Example:

```r
processed <- Cardinal::peakProcess(
    msData,
    method = "diff",
    SNR = 3,
    sampleSize = 0.1,
    filterFreq = 0.01,
    tolerance = 0.5,
    units = "mz"
)
```

This can combine:

- Peak detection.
- Reference peak identification.
- Peak alignment.
- Peak-frequency filtering.

When `sampleSize` is supplied, reference peaks can be estimated from a subset of spectra.

The reference peaks are then used to summarize the full dataset.

### Implications for moleculaR

This is potentially useful for generating a compact representation for repeated large-scale analyses.

However, the preprocessing may remove rare analytes.

Therefore, it should remain optional and must not replace direct centroided-spectrum extraction by default.

## 8.6 `bin()`

Cardinal provides:

```r
Cardinal::bin()
```

Example:

```r
binned <- Cardinal::bin(
    msData,
    resolution = 10,
    units = "ppm",
    method = "sum"
)
```

Supported aggregation methods include:

- `sum`
- `mean`
- `max`
- `min`
- Interpolation methods.

The resulting dataset has a shared reference axis.

### Important distinction

Binning is not equivalent to searching for a theoretical analyte using an FWHM-dependent Gaussian.

Cardinal's `bin()` should not be used as a substitute for moleculaR's Gaussian intensity extraction.

Use it when an aligned or resampled dataset is needed.

## 8.7 `spectrapply()`

Cardinal provides:

```r
Cardinal::spectrapply()
```

This is particularly important for moleculaR.

It allows a function to be applied to individual spectra while Cardinal manages iteration and chunked access.

Example:

```r
result <- Cardinal::spectrapply(
    msData,
    function(intensity, mz) {
        sum(intensity[mz > 700 & mz < 800])
    },
    index = "mz",
    simplify = TRUE
)
```

The first argument to the supplied function receives the spectral intensities.

When `index = "mz"` is provided, the second argument receives the corresponding m/z values.

### Why this is important

It provides a potential implementation route for analyte extraction without creating a global sparse matrix.

The agent must investigate its behavior on both:

- `MSImagingArrays`.
- `MSImagingExperiment`.

For the initial implementation, use serial execution.

Later, benchmark parallel and custom chunk-based implementations.

## 8.8 Chunking

Cardinal and matter provide chunk-based operations.

The following Cardinal options are relevant:

```r
Cardinal::getCardinalNChunks()
```

```r
Cardinal::setCardinalNChunks()
```

```r
Cardinal::getCardinalChunksize()
```

```r
Cardinal::setCardinalChunksize()
```

Many Cardinal functions also accept:

```r
chunkopts = list()
```

The exact chunk configuration should be determined experimentally.

### Important requirement

Chunk sizes must not be hard-coded without justification.

Performance depends on:

- Number of pixels.
- Peaks per spectrum.
- Storage device.
- Available RAM.
- Number of target analytes.
- Parallel processing configuration.

Benchmark different chunk sizes.

## 8.9 Parallel processing

Cardinal supports BiocParallel.

Recommended initial configuration:

```r
BPPARAM <- BiocParallel::SerialParam()
```

For supported parallel execution:

```r
BPPARAM <- BiocParallel::MulticoreParam(
    workers = 4
)
```

or:

```r
BPPARAM <- BiocParallel::SnowParam(
    workers = 4
)
```

The latter is useful for cross-platform testing.

Do not assume Unix-style forked parallelization is available on Windows.

### File-backed data considerations

Parallel workers must have access to the underlying spectral files.

Investigate:

- File path accessibility.
- Serialization overhead.
- Number of simultaneous file readers.
- Memory duplication.
- Disk contention.
- Error recovery.
- Reproducibility.

Parallelization should only be enabled when benchmarks demonstrate a benefit.

---

# 9. Function-by-Function Migration Map

The following table defines the intended migration.

| Existing function | Intended Cardinal implementation | Action |
|---|---|---|
| `readCentrData()` | `readMSIData()` | Replace or deprecate |
| `readSingleSpect()` | Numeric spectrum representation | Rewrite |
| `createSparseMat()` | Native Cardinal data access | Remove from main workflow |
| `searchAnalyte()` | Direct Cardinal intensity extraction | Major rewrite |
| `batchAnalyteSearch()` | Chunked multi-analyte extraction | Major optimization |
| `batchLipidSearch()` | New batch extraction engine | Refactor |
| `estimateFwhm()` | Cardinal-compatible spectra and numeric processing | Rewrite MALDIquant-dependent parts |
| `useDefaultFwhm()` | Existing models | Preserve |
| `getFwhm()` | Existing implementation | Preserve |
| `gaussWeight()` | Existing implementation | Preserve |
| `normalizeIntensity()` | Cardinal normalization + custom methods | Refactor |
| `filterPeaks()` | Cardinal feature filtering | Refactor |
| `createSpatialWindow()` | Cardinal coordinates | Minor adaptation |
| `analytePointPattern()` | Existing S3 class | Preserve |
| `ionImage2spp()` | Existing function | Preserve |
| `probMap()` | Existing spatial methods | Preserve |
| `superimposeAnalytes()` | Existing spatial methods | Preserve |
| `subsetAnalytes()` | Existing spatial methods | Preserve |
| `recreateExData()` | Cardinal-based example-data workflow | Rewrite |
| Shiny applications | Direct Cardinal objects | Refactor |

Each migration requires its own tests.

The following sections specify the most important implementations.

---

# 10. Replacing `readCentrData()`

## 10.1 Current behavior

`readCentrData()` currently:

1. Reads the beginning of the imzML XML file.
2. Checks whether the file appears to use continuous representation.
3. Rejects continuous-format input.
4. Calls `MALDIquantForeign::importImzMl()`.
5. Returns a list of MALDIquant `MassPeaks` objects.

This approach should be replaced.

## 10.2 Proposed public function

Introduce:

```r
readMsi()
```

Suggested interface:

```r
readMsi <- function(
    path,
    memory = FALSE,
    as = "MSImagingArrays",
    ...
) {
    Cardinal::readMSIData(
        file = path,
        memory = memory,
        as = as,
        ...
    )
}
```

This is a starting point, not the complete implementation.

Additional validation should check:

- File existence.
- imzML/ibd pairing.
- Supported formats.
- Coordinate availability.
- Spectrum representation.
- Dataset integrity.
- Whether the dataset contains multiple runs.

## 10.3 Compatibility

The old `readCentrData()` may be retained temporarily as a deprecated wrapper returning a Cardinal object.

However, because its return type changes substantially, the migration is not fully backward-compatible.

The documentation must explain this explicitly.

## 10.4 Required tests

Test:

- Continuous imzML import.
- Processed imzML import.
- Centroided imzML.
- Processed but non-centroided imzML.
- Missing imzML file.
- Missing ibd file.
- Corrupt binary data.
- Invalid coordinates.
- Repeated coordinates.
- Multiple experimental runs.
- File-backed import.
- In-memory import.
- Metadata preservation.

Use `CardinalIO::exampleImzMLFile()` for basic import tests where appropriate.

---

# 11. Replacing `createSparseMat()`

## 11.1 Design decision

Do not replace `createSparseMat()` with a function that immediately constructs another complete matrix.

The primary Cardinal workflow should not require an equivalent global sparse matrix.

Instead, analyte extraction should operate directly on Cardinal objects.

## 11.2 Legacy function

The old `createSparseMat()` may remain during migration for reference tests.

It should eventually be removed from the main workflow.

If retained as a compatibility utility, it must not be invoked automatically.

## 11.3 Required tests

Verify that:

- A Cardinal object can be passed directly to `searchAnalyte()`.
- The complete spectral dataset is not converted into an ordinary matrix.
- The number of spectra is preserved.
- Pixel coordinates remain correctly associated with spectra.
- Analyte extraction works without an aligned global mass axis.

The absence of a global sparse matrix conversion is an explicit architectural acceptance criterion.

---

# 12. Reimplementing `searchAnalyte()`

## 12.1 Priority

**Highest priority.**

This is the central scientific and technical task of the migration.

The success of the entire project depends on this implementation.

## 12.2 Required API

Preserve the existing public interface where feasible:

```r
searchAnalyte(
    m,
    fwhm,
    spData,
    wMethod = "Gaussian",
    spwin = NA,
    verifiedMasses = NA,
    confirmedOnly = FALSE,
    metaData = list()
)
```

The `spData` argument should accept Cardinal objects.

The argument name may be retained initially to minimize disruption.

A future API revision could introduce a more intuitive `msData` argument.

Do not change the argument name and its semantics simultaneously without a compatibility plan.

## 12.3 Current scientific logic

For a theoretical analyte mass \(m\), moleculaR uses the estimated FWHM to determine a Gaussian standard deviation:

\[
\sigma = \frac{\mathrm{FWHM}(m)}{2\sqrt{2\ln 2}}
\]

The nominal search tolerance is:

\[
\Delta m = 3\sigma
\]

The Gaussian weight is:

\[
w(m_j)=
\exp\left(
-\frac{(m_j-m)^2}{2\sigma^2}
\right)
\]

This is equivalent to the normalized Gaussian weighting implemented by `gaussWeight()`.

For the Gaussian aggregation method:

\[
I_p(m)=\sum_j I_{pj}w(m_{pj})
\]

where the sum is taken over the selected measured peaks.

## 12.4 Important legacy behavior

The existing implementation does not simply select every mass satisfying:

```r
abs(mz - m) <= 3 * sigma
```

Instead, it uses:

```r
MALDIquant::match.closest()
```

to identify:

- The nearest feature to the target mass.
- The lower feature index.
- The upper feature index.

It subsequently extracts the complete range of feature columns between those indices.

The resulting interval may differ from a strict mathematical mass window.

This distinction is particularly important near search boundaries.

### Required action

Document and test the difference between:

1. Legacy nearest-index selection.
2. Strict mass-window selection.

For the new Cardinal implementation, strict FWHM-based window selection is the preferred default scientific definition.

However, do not claim numerical equivalence with the legacy implementation until the difference has been measured.

A short-lived internal compatibility implementation may be useful during validation.

It does not need to become a permanent user-facing backend option.

## 12.5 Gaussian aggregation

The default method should calculate:

```r
sum(intensity[selected] * weights)
```

where:

```r
weights <- gaussWeight(
    x = mz[selected],
    m = m,
    fwhm = fwhm
)
```

The selected peaks must correspond to the same spectrum.

The calculated intensities must remain in the original pixel order.

## 12.6 Alternative aggregation methods

Preserve:

```r
"Gaussian"
"sum"
"mean"
"max"
```

### Sum

```r
sum(intensity[selected])
```

### Maximum

```r
max(intensity[selected])
```

Handle empty selections explicitly.

### Mean

The mean method requires special attention.

The existing implementation calculates a mean over columns in the global sparse intensity matrix.

Missing peaks are represented by zeros.

Therefore, the denominator is the number of selected global feature columns, not necessarily the number of observed peaks within a particular pixel.

A direct per-spectrum implementation using:

```r
mean(intensity[selected])
```

can produce different results.

### Required decision

The agent must explicitly document whether the new mean method:

- Preserves the global-feature denominator.
- Uses only observed peaks.
- Uses a defined fixed reference grid.

The preferred architecture should avoid reconstructing a global mass axis merely to support this method.

However, changing the meaning of `mean` must be treated as a documented scientific/API change.

Do not silently replace the old denominator.

## 12.7 Initial implementation using `spectrapply()`

A conceptual starting point is:

```r
.extractIonArrays <- function(
    x,
    m,
    fwhm,
    method = "Gaussian"
) {

    sigma <- fwhm / (2 * sqrt(2 * log(2)))

    tolerance <- 3 * sigma

    Cardinal::spectrapply(
        x,
        FUN = function(intensity, mz) {

            selected <- which(
                abs(mz - m) <= tolerance
            )

            if (length(selected) == 0L) {
                return(0)
            }

            selectedMz <- mz[selected]

            selectedIntensity <- intensity[selected]

            switch(
                method,

                "Gaussian" = {
                    weights <- gaussWeight(
                        x = selectedMz,
                        m = m,
                        fwhm = fwhm
                    )

                    sum(selectedIntensity * weights)
                },

                "sum" = sum(selectedIntensity),

                "mean" = mean(selectedIntensity),

                "max" = max(selectedIntensity),

                stop("Unsupported aggregation method.")
            )
        },
        index = "mz",
        simplify = TRUE,
        BPPARAM = BiocParallel::SerialParam()
    )
}
```

**This code is a prototype, not a production implementation.**

In particular:

- The mean method does not reproduce the legacy denominator.
- Input validation is incomplete.
- Missing and non-finite intensity handling is incomplete.
- Processing-state handling is not implemented.
- Performance has not been optimized.
- The implementation has not been executed against the selected Cardinal release.

The agent must address these issues before integrating it into the package.

## 12.8 Supporting `MSImagingExperiment`

Implement extraction for shared-axis data.

The relevant m/z values can be obtained using:

```r
Cardinal::mz(x)
```

Select the corresponding feature rows.

Then aggregate their intensities across pixels.

Avoid converting the complete intensity matrix to an ordinary R matrix.

For example, only the relevant feature rows should be retrieved.

A small illustrative expression is:

```r
mzAxis <- Cardinal::mz(x)

selected <- which(
    abs(mzAxis - m) <= tolerance
)

selectedData <- Cardinal::intensity(x)[
    selected,
    ,
    drop = FALSE
]
```

The aggregation must then be performed using a method compatible with the actual matrix representation.

The implementation should use chunked or optimized matrix operations where beneficial.

## 12.9 Returning point patterns

The extracted intensities must be combined with:

```r
Cardinal::coord(x)
```

The current `searchAnalyte()` includes only pixels satisfying:

```r
combinedIntensity > 0
```

Preserve this behavior.

The returned object must remain compatible with:

```r
analytePointPattern()
```

and:

```r
probMap()
```

Preserve:

- x coordinates.
- y coordinates.
- Intensity marks.
- Analyte m/z.
- Verification metadata.
- User-supplied metadata.
- Spatial window.
- Expected class structure.

## 12.10 Handling zero intensities

Distinguish between:

1. A measured pixel with zero analyte intensity.
2. A pixel absent from the MSI acquisition.
3. A missing or invalid intensity value.

These are not interchangeable.

Zero-intensity acquired pixels must remain part of the tissue window.

They should be omitted from the detected analyte points, consistent with the existing default.

Missing pixels should not be fabricated.

## 12.11 Handling multiple runs

The current implementation assumes a single spatial analysis domain.

Cardinal can contain several experimental runs.

Different runs may contain identical x and y coordinates.

Do not combine these runs into one point pattern without an explicit selection or coordinate transformation.

### Initial recommended behavior

If a Cardinal object contains multiple runs, require the caller to select a run or subset before generating a single MPM.

A future API could support an explicit `run` argument.

The first implementation should prioritize correctness over convenience.

## 12.12 Verification against METASPACE

The existing implementation uses:

```r
MALDIquant::match.closest()
```

to match theoretical masses to externally verified masses.

Replace this with an independent numerical matching function.

The replacement must:

- Accept numeric verified masses.
- Apply the correct absolute tolerance.
- Select the nearest verified mass.
- Handle empty input.
- Handle ties deterministically.
- Handle invalid values.
- Preserve the meaning of `mzConfirmed`.
- Preserve `confirmedOnly`.

Do not change how annotations are interpreted during the infrastructure migration.

## 12.13 Required tests

Construct synthetic spectra containing:

- One exact matching peak.
- One peak within tolerance.
- One peak outside tolerance.
- Multiple matching peaks.
- Peaks at the tolerance boundary.
- No matching peaks.
- Zero-intensity peaks.
- Closely spaced peaks.
- Duplicate m/z values.
- Empty spectra.
- Mass values near the global mass-range boundary.

For each spectrum, calculate the expected intensity manually.

Test all four aggregation methods.

Verify that the final intensity vector is associated with the correct coordinates.

---

# 13. Reimplementing `batchAnalyteSearch()`

## 13.1 Current behavior

The current function repeatedly calls:

```r
searchAnalyte()
```

using:

```r
parallel::mclapply()
```

This is appropriate when a global sparse matrix is already available.

It may become inefficient for file-backed data.

Repeatedly scanning a large imzML dataset for thousands of analytes can create excessive disk I/O.

## 13.2 Target behavior

The new implementation should evaluate multiple target masses during a single traversal of the spectral data whenever practical.

Conceptual workflow:

```text
Target masses
      |
      v
Sort and validate targets
      |
      v
Determine FWHM for each target
      |
      v
Construct mass search intervals
      |
      v
Read a chunk of spectra
      |
      v
Evaluate all relevant targets
      |
      v
Store extracted intensities
      |
      v
Continue with next chunk
      |
      v
Combine results
      |
      v
Construct analyte point patterns
```

## 13.3 Suggested internal function

Implement:

```r
.msiExtractIons(
    x,
    mz,
    fwhm,
    method = "Gaussian",
    chunkSize = NULL,
    BPPARAM = BiocParallel::SerialParam()
)
```

The function should return a compact intensity representation.

Conceptually:

```text
                   Analyte 1   Analyte 2   Analyte 3

Pixel 1              10.2        0.0        15.3
Pixel 2               0.0       12.4         0.0
Pixel 3               8.1        4.2         0.0
```

Choose the actual storage representation based on memory requirements.

Do not assume a dense pixel-by-analyte matrix will always be appropriate.

## 13.4 Optimization strategies

Investigate:

- Sorting target masses.
- Binary-search-based peak lookup.
- Reusing mass-range calculations.
- Vectorized Gaussian weighting.
- Chunked processing.
- Sparse result storage.
- Reusing extracted ion intensities.
- Grouping overlapping mass windows.
- Efficient allocation of output arrays.

Avoid repeated memory allocation inside deeply nested loops.

## 13.5 Caching

Introduce optional caching of extracted analyte intensities.

The cache key should include:

- Dataset identity.
- Processing state.
- Target masses.
- FWHM model.
- Aggregation method.
- Mass-selection rule.
- Run or ROI selection.

A cache must not return intensities calculated from a different normalization state or dataset.

## 13.6 Required tests

Verify that:

```r
batchAnalyteSearch()
```

produces the same per-analyte results as repeated calls to:

```r
searchAnalyte()
```

using the same Cardinal object.

Test:

- Duplicate target masses.
- Unsorted target masses.
- Missing masses.
- Out-of-range masses.
- Empty target vectors.
- Annotation filtering.
- Single-core execution.
- Parallel execution.
- Deterministic output ordering.

Performance should be evaluated only after numerical correctness is established.

---

# 14. Reimplementing `batchLipidSearch()`

## 14.1 Existing behavior

The current implementation searches the SwissLipids database and evaluates several possible adducts.

Supported adducts include:

```r
"M-H"
"M+H"
"M+Na"
"M+K"
```

The implementation generates metadata describing the detected lipid species.

It subsequently combines the detected analytes into the appropriate spatial structures.

## 14.2 Target design

Separate lipid database processing from spectral extraction.

Recommended architecture:

```text
SwissLipids database
         |
         v
Select lipid species
         |
         v
Generate theoretical adduct masses
         |
         v
Apply annotation filtering
         |
         v
Batch spectral extraction
         |
         v
Attach lipid metadata
         |
         v
Combine analyte point patterns
```

Do not implement a separate spectral extraction algorithm for each adduct.

All theoretical masses should be passed through the same optimized extraction infrastructure.

## 14.3 Preserve metadata

The current code attaches information such as:

- Lipid identifier.
- Ionization mode.
- Adduct.
- Sum formula.
- Abbreviation.
- Number of double bonds.
- Lipid class.
- Chain length.

Preserve these fields and their relationships with the returned point patterns.

## 14.4 Required tests

Test every supported adduct.

Verify:

- Correct theoretical masses.
- Correct annotation matching.
- Correct lipid metadata.
- Correct positive/negative ionization assignment.
- Correct exclusion of missing masses.
- Correct handling of repeated masses.
- Correct ordering and grouping.
- Equivalence to the new single-analyte extraction.
- Compatibility with collective projection functions.

Database processing should not require loading the MSI spectral data repeatedly.

---

# 15. Reimplementing FWHM Estimation

## 15.1 Importance

FWHM estimation is central to the scientific methodology of moleculaR.

The FWHM determines the mass tolerance used during analyte detection.

Any numerical change in FWHM can alter:

- Which peaks are included.
- Gaussian weights.
- Extracted intensities.
- Detected pixels.
- Molecular probability maps.

Therefore, this functionality requires especially careful validation.

## 15.2 Functions to preserve

Preserve the public functions:

```r
estimateFwhm()
```

```r
useDefaultFwhm()
```

```r
getFwhm()
```

Preserve the existing S3 `fwhm` object structure wherever possible.

## 15.3 Existing input types

The current `estimateFwhm()` supports:

- A single MALDIquant `MassSpectrum`.
- A list of MALDIquant `MassSpectrum` objects.
- A list of `MassPeaks` objects with stored FWHM metadata.
- An instrument-dependent default model when no spectra are supplied.

These inputs are handled through MALDIquant-specific code.

## 15.4 Proposed new input types

Support:

1. A Cardinal MSI object containing profile spectra.
2. A numeric m/z vector and an intensity vector.
3. A list or data frame containing m/z and intensity arrays.
4. Precomputed peak positions and widths.
5. Existing instrument-based models.

A single profile spectrum imported from TSV should not require a MALDIquant `MassSpectrum` object.

## 15.5 Reimplementing width estimation

The existing implementation:

1. Detects peaks.
2. Determines the relevant peak indices.
3. Calculates half-height values.
4. Searches to the left and right of each peak.
5. Estimates peak widths.
6. Fits an FWHM-versus-m/z relationship.

The numerical width-estimation algorithm can be adapted to ordinary numeric vectors.

Conceptual internal interface:

```r
.estimatePeakWidths(
    mz,
    intensity,
    peakIndices
)
```

The agent must inspect the current `.fwhm()` and `.getFwhm()` implementations carefully.

Do not replace the existing numerical method with a different width estimator without validation.

## 15.6 Peak detection

Cardinal provides peak-picking functionality.

However, Cardinal's peak-picking algorithms are not automatically equivalent to MALDIquant's `detectPeaks(method="SuperSmoother")`.

Differences in peak detection can change the fitted FWHM curve.

### Recommended approach

Preserve the existing FWHM fitting methods first.

Replace the input handling and peak-width calculation.

Then investigate whether Cardinal's peak detection can serve as a suitable replacement.

Document any resulting differences.

## 15.7 Default FWHM models

The existing `useDefaultFwhm()` implements instrument-dependent models.

These include:

- FTICR.
- TOF-reflector.
- TIMS-TOF.
- User-specified resolving power.

These models should be preserved.

### Important issue: bundled FTICR reference data

The current `useDefaultFwhm()` loads the package's bundled `processed-example-Data` dataset to obtain an existing reference `fwhmObj`.

This reference is used when constructing the default nonlinear FTICR model.

Therefore, removing the legacy example dataset or modifying its contents could unintentionally change the default FWHM behavior.

### Required action

Preserve the numerical FWHM reference independently of MALDIquant.

A suitable approach is to store the reference m/z values, FWHM values, and required interpolation information as native R data.

The default model must reproduce the existing results at representative reference masses.

## 15.8 Required tests

Test:

- Constant-width model.
- FTICR model.
- TIMS-TOF model.
- TOF-reflector model.
- User-specified resolving power.
- Interpolation.
- Extrapolation.
- Invalid resolution.
- Invalid spectrum.
- Missing peaks.
- Overlapping peaks.
- Insufficient peak counts.
- Profile spectra with known synthetic peak widths.

Compare the new implementation with the legacy results.

Where the peak-detection method differs, report both the differences in selected peaks and the resulting FWHM estimates.

---

# 16. Reimplementing `readSingleSpect()`

## 16.1 Existing behavior

The current function imports a two-column TSV file and creates a MALDIquant `MassSpectrum`.

This is primarily used for FWHM estimation.

## 16.2 Proposed replacement

Return a lightweight numeric representation.

For example:

```r
readSingleSpect <- function(path) {

    df <- utils::read.delim(path)

    spectrum <- list(
        mz = as.numeric(df[[1]]),
        intensity = as.numeric(df[[2]])
    )

    class(spectrum) <- "moleculaRSpectrum"

    spectrum
}
```

This is a suggested starting point.

Add appropriate input validation.

The new object should not depend on MALDIquant.

## 16.3 Required tests

Check:

- Two-column TSV input.
- Correct numeric conversion.
- Valid m/z ordering.
- Missing values.
- Duplicate m/z values.
- Empty spectra.
- Incorrect number of columns.
- Compatibility with `estimateFwhm()`.

Document the changed return type.

---

# 17. Reimplementing Normalization

## 17.1 Existing methods

moleculaR currently supports:

```r
"medianFoldChange"
"RMS"
"TIC"
"median"
"IS"
```

These must remain available.

## 17.2 TIC normalization

The existing implementation divides each intensity by its spectrum's total intensity.

Cardinal provides:

```r
Cardinal::normalize(
    msData,
    method = "tic"
)
```

Compare the exact resulting intensities.

Verify the role of Cardinal's `scale` argument.

## 17.3 RMS normalization

The existing RMS normalization uses:

\[
\mathrm{RMS}=\sqrt{\frac{1}{n}\sum_i I_i^2}
\]

and divides the spectral intensities by this factor.

Cardinal provides:

```r
Cardinal::normalize(
    msData,
    method = "rms"
)
```

Verify numerical equivalence.

## 17.4 Median normalization

The existing moleculaR method divides intensities by the median spectral intensity.

Implement this method for Cardinal objects if no equivalent method is available with the required semantics.

Use chunked or deferred processing.

## 17.5 Internal-standard normalization

The current implementation:

1. Finds a specified internal-standard peak.
2. Retrieves its intensity in each spectrum.
3. Replaces missing internal-standard intensities with the median of available intensities.
4. Divides the spectrum by the resulting normalization factor.

This behavior must be preserved unless there is a documented reason to change it.

Cardinal's `reference` normalization may provide similar functionality, but it must not be assumed to reproduce the same missing-value policy.

Implement a custom method when necessary.

## 17.6 Median fold-change normalization

The current implementation uses `MALDIquant::mergeMassPeaks()` to construct a reference spectrum.

It then calculates spectrum-specific median fold changes.

This method depends on shared feature positions.

For Cardinal, a shared-axis representation is likely to be appropriate.

The implementation must specify:

- How the reference spectrum is constructed.
- How features are matched.
- How missing peaks are handled.
- How zero reference intensities are handled.
- How normalization factors are calculated.

Preserve the existing numerical definition where feasible.

## 17.7 Custom deferred preprocessing

Cardinal's `addProcessing()` can be investigated for implementing custom normalization steps.

It accepts a function operating on spectral intensities and corresponding locations.

For example, a custom median normalization may be implemented without manually transforming the entire dataset into an R list.

The exact implementation must be validated against the installed Cardinal release.

## 17.8 Required tests

For each normalization method:

- Compare normalization factors.
- Compare normalized intensities.
- Test zero-intensity spectra.
- Test missing reference peaks.
- Test non-finite values.
- Test different spectrum lengths.
- Test file-backed processing.
- Test repeated normalization.
- Test deferred processing state.

Normalization must not silently introduce `NaN` or `Inf` values.

---

# 18. Reimplementing `filterPeaks()`

## 18.1 Existing behavior

The current function removes peaks occurring in fewer than a specified proportion of pixels.

It uses the sparse intensity matrix to calculate feature frequencies.

The filtering is then applied to the original MALDIquant spectra.

## 18.2 Cardinal replacement

For aligned Cardinal data, peak frequencies can be represented in:

```r
Cardinal::featureData(msData)
```

Following peak alignment, Cardinal may provide:

```r
count
freq
```

Filtering can be performed using:

```r
filtered <- Cardinal::subsetFeatures(
    msData,
    freq >= minFreq
)
```

The exact threshold semantics must match the intended moleculaR behavior.

## 18.3 Important limitation

For unaligned `MSImagingArrays`, there is no single shared feature axis.

Therefore, global peak-frequency filtering requires defining how peaks are matched across spectra.

This necessarily introduces an alignment or mass-tolerance decision.

### Recommended behavior

Apply global peak-frequency filtering primarily to aligned data.

Do not require this step for direct single-analyte extraction.

## 18.4 Required tests

Verify:

- Frequency calculations.
- Threshold behavior.
- Peaks exactly at the threshold.
- Peaks detected in only one spectrum.
- Duplicate peak assignments.
- Missing intensities.
- Preservation of metadata.
- Correct behavior when all peaks are removed.

---

# 19. Spatial Data and Coordinate Handling

## 19.1 Cardinal coordinates

Use:

```r
Cardinal::coord(msData)
```

Convert to an ordinary two-column data frame when passing the information to the existing spatial functions.

Example:

```r
coordinates <- as.data.frame(
    Cardinal::coord(msData)
)

coordinates <- coordinates[
    ,
    c("x", "y"),
    drop = FALSE
]
```

Verify that the output includes the expected coordinate names.

## 19.2 Tissue window construction

Preserve:

```r
createSpatialWindow()
```

This function should receive the coordinates of all acquired tissue pixels.

It must not receive only the pixels in which a particular analyte was detected.

Otherwise, the spatial reference domain would change between analytes.

## 19.3 Coordinate orientation

Cardinal and spatstat may use different conventions for presenting spatial data.

The migration must preserve:

- Pixel ordering.
- x/y orientation.
- Tissue shape.
- Pixel spacing.
- Missing pixels.
- Disconnected tissue regions.
- Spatial boundaries.

Test asymmetric tissue shapes to detect accidental transposition or mirroring.

## 19.4 Physical pixel size

Where available, preserve physical pixel-size metadata.

Distinguish between:

- Pixel indices.
- Physical distances.
- Spatial analysis units.

Do not silently apply physical scaling to the existing MPM algorithms.

Any introduction of physical spatial units should be implemented as a separate documented feature.

## 19.5 Multiple runs

Multiple experimental runs must be treated as separate spatial domains unless explicitly combined using a valid coordinate transformation.

Never combine runs solely because their x/y coordinates have the same names.

---

# 20. Preserving `analytePointPattern`

The existing `analytePointPattern` class is the interface between the MSI processing layer and the spatial analysis layer.

Its structure must remain compatible with existing methods.

Important components include:

- spatstat `ppp` representation.
- Coordinates.
- Spatial window.
- Intensity marks.
- Analyte identifiers.
- Associated metadata.

The current point-pattern constructor also generates an internal identifier linking point marks to analyte metadata.

The agent must preserve this relationship.

### Required validation

For the same input intensity vector and coordinates, the Cardinal implementation should construct the same logical point pattern as the legacy implementation.

Compare the meaningful scientific fields rather than incidental identifiers generated from time or randomness.

Do not compare entire serialized objects when nondeterministic metadata are present.

---

# 21. Preserving Molecular Probability Mapping

## 21.1 `probMap()`

The primary spatial probability mapping function should remain unchanged during the initial migration.

It accepts `analytePointPattern` data rather than MALDIquant spectra.

Therefore, the correct migration boundary is before this function.

## 21.2 Scientific validation

Compare MPM outputs generated from equivalent analyte point patterns.

Important components include:

- Bandwidth.
- Molecular density estimate.
- CSR reference representation.
- Reference density estimate.
- Statistical significance.
- Hotspot locations.
- Coldspot locations.
- Hotspot boundaries.
- Coldspot boundaries.

Use fixed random seeds.

Use deterministic reference inputs where possible.

## 21.3 Avoid scope expansion

Do not simultaneously redesign:

- Kernel density estimation.
- CSR simulation.
- Spatial autocorrelation.
- Statistical testing.
- FDR correction.

These should remain independent of the MSI infrastructure migration.

If performance issues are discovered in these methods, document them for a separate optimization stage.

---

# 22. Shiny Application Migration

## 22.1 Existing assumptions

The current Shiny application assumes that initialization generates:

- A list of MALDIquant spectra.
- A global sparse intensity matrix.
- A shared `mzAxis`.
- A corresponding coordinate representation.

The application directly invokes MALDIquant processing methods.

These assumptions must be removed.

## 22.2 Proposed server-side state

Maintain:

```r
msData
```

as the primary Cardinal object.

Additional state may include:

```r
fwhmObj
spatialWindow
selectedRun
processingSettings
ionCache
```

Do not maintain a duplicate complete spectral representation.

## 22.3 Updated initialization workflow

```text
Upload imzML and ibd
        |
        v
Validate files
        |
        v
Cardinal import
        |
        v
Validate spectrum representation
        |
        v
Select run / tissue region
        |
        v
Optional preprocessing
        |
        v
Create spatial window
        |
        v
Initialize FWHM model
        |
        v
Ready for MPM analysis
```

## 22.4 Mass selection

The old Shiny application selects a nearest mass using:

```r
spData$mzAxis
```

This is inappropriate for unaligned `MSImagingArrays`.

The new application should allow the user to enter a theoretical mass directly.

The mass must be validated against the dataset's supported mass range.

The extraction method should then evaluate the required FWHM-based window.

Do not require a shared feature axis merely to validate the user-entered mass.

## 22.5 File-backed storage lifecycle

Cardinal objects may reference the original imzML/ibd files.

These files must remain accessible while the data are being analyzed.

For Shiny uploads:

- Copy uploaded files to a controlled location.
- Preserve the imzML/ibd relationship.
- Do not rely on temporary upload paths after their lifecycle ends.
- Use session-specific storage where appropriate.
- Avoid filename collisions.
- Clean up files after they are no longer required.
- Ensure users cannot access other sessions' uploaded data.

A serialized Cardinal object is not necessarily self-contained.

## 22.6 User feedback

Large MSI processing operations may take significant time.

Provide progress indicators for:

- File validation.
- Dataset import.
- Preprocessing.
- FWHM estimation.
- Batch extraction.
- MPM calculation.

Do not freeze the interface unnecessarily.

## 22.7 Caching

Cache extracted ion intensities when beneficial.

Invalidate the cache whenever relevant processing settings change.

Do not cache results solely by theoretical m/z.

## 22.8 Required Shiny tests

Verify:

- Successful file upload.
- Missing ibd handling.
- Invalid input.
- Centroided-data handling.
- Profile-data rejection or preprocessing.
- FWHM initialization.
- Single-analyte MPM generation.
- Batch analyte analysis.
- Lipid search.
- Visualization.
- Downloaded results.
- File cleanup.
- Multiple sessions.
- Large dataset memory behavior.

Where practical, use `shinytest2` for automated interface tests.

---

# 23. Example Data and Vignette Migration

## 23.1 Current example data

The existing example-data workflow is implemented through:

```r
recreateExData()
```

It relies on MALDIquant import, binning, filtering, and sparse matrix construction.

This must be updated.

## 23.2 New example-data strategy

Use a small Cardinal-compatible example dataset.

The dataset should be suitable for:

- Installation tests.
- Vignette execution.
- Single-analyte extraction.
- FWHM estimation.
- Batch searches.
- MPM generation.

Avoid packaging unnecessarily large binary files.

A small in-memory Cardinal object is suitable for unit tests.

A small imzML/ibd pair is useful for import and file-backed tests.

## 23.3 Legacy data considerations

Do not immediately delete the existing example data.

The bundled `fwhmObj` is currently used by the default FTICR model.

Preserve the required FWHM information first.

## 23.4 Vignettes

Update the package walkthrough to demonstrate:

```r
library(moleculaR)
library(Cardinal)
```

```r
msData <- readMsi(
    "dataset.imzML",
    memory = FALSE
)
```

```r
fwhmObj <- useDefaultFwhm(
    modality = "FTICR"
)
```

```r
analyte <- searchAnalyte(
    m = 788.5447,
    fwhm = getFwhm(fwhmObj, 788.5447),
    spData = msData
)
```

```r
mpm <- probMap(
    analyte,
    seed = 1
)
```

These examples represent the intended API and must be tested after implementation.

Provide separate examples for:

- Direct centroided-spectrum extraction.
- Optional Cardinal preprocessing.
- Batch analyte analysis.
- SwissLipids searches.
- File-backed data handling.

---

# 24. Dependency and Package Infrastructure Changes

## 24.1 `DESCRIPTION`

The existing package imports MALDIquant.

The final package should declare Cardinal as a dependency.

The precise dependency arrangement must be determined after implementation.

Expected changes include:

- Adding Cardinal.
- Adding required Bioconductor infrastructure dependencies where directly used.
- Removing MALDIquant.
- Removing MALDIquantForeign.
- Updating the minimum supported R version.
- Updating package documentation.
- Updating suggested testing dependencies.

Do not add matter or BiocParallel to direct imports unless moleculaR actually calls their APIs.

Avoid unnecessary dependency duplication.

## 24.2 Bioconductor version compatibility

Cardinal releases are tied to compatible Bioconductor and R versions.

Pin a supported environment.

Do not mix arbitrary package versions from different Bioconductor releases.

Use:

```r
BiocManager::valid()
```

to inspect the installed environment.

Record:

```r
sessionInfo()
```

in the validation reports.

## 24.3 `renv`

The repository uses `renv`.

Update the dependency lockfile only after verifying compatibility.

Do not restore an old dependency environment that prevents the selected Cardinal release from installing.

The migration must produce a reproducible installation process.

## 24.4 `NAMESPACE`

Update exports and imports using the project's roxygen2 workflow.

Do not manually modify a generated `NAMESPACE` without regeneration.

## 24.5 `Collate`

The current `DESCRIPTION` contains an explicit `Collate` list.

Any newly created R source files must be integrated correctly.

Otherwise, the package may fail during installation or load functions in the wrong order.

Either maintain the list correctly or adopt a deliberate alternative that preserves source ordering.

## 24.6 Package checks

Run:

```r
devtools::test()
```

```r
devtools::check()
```

and:

```r
BiocManager::valid()
```

The final implementation should pass R package checks without requiring MALDIquant.

---

# 25. Scientific Validation Strategy

This is a mandatory component of the migration.

The development agent must not rely only on successful function execution.

The results must be scientifically validated.

## 25.1 Validation levels

Use four validation levels.

### Level 1: Input data

Compare:

- Spectrum count.
- m/z values.
- Intensity values.
- Coordinates.
- Experimental runs.
- Centroid status.
- Metadata.

### Level 2: Analyte extraction

Compare:

- Selected peaks.
- Gaussian weights.
- Aggregated intensities.
- Detected pixel counts.
- Detected pixel coordinates.
- Analyte metadata.

### Level 3: Spatial analysis

Compare:

- Point patterns.
- KDE bandwidths.
- Density images.
- CSR reference outputs.
- Hotspot/coldspot classifications.
- Final probability maps.

### Level 4: Computational performance

Compare:

- Import time.
- Peak memory usage.
- Preprocessing time.
- Single-analyte extraction time.
- Batch extraction time.
- MPM calculation time.
- Disk I/O where measurable.

## 25.2 Establishing a legacy reference

Before modifying the existing functions, run the legacy implementation on a fixed reference dataset.

Save the scientifically relevant outputs.

For example:

```text
tests/fixtures/legacy/
    dataset-metadata.rds
    reference-masses.rds
    reference-fwhm.rds
    extracted-intensities.rds
    analyte-coordinates.rds
    point-pattern-summary.rds
    mpm-summary.rds
```

Avoid saving enormous or unnecessary serialized objects.

Prefer compact numerical reference data.

## 25.3 Numerical tolerances

Use strict numerical comparisons where the mathematical operations are intended to be identical.

For example:

```r
testthat::expect_equal(
    actual,
    expected,
    tolerance = 1e-10
)
```

However, distinguish between:

- Floating-point differences.
- Different peak-selection rules.
- Different preprocessing algorithms.
- Scientifically meaningful differences.

A universal tolerance is not appropriate for every result.

Use absolute and relative errors where necessary.

Record why each tolerance was chosen.

## 25.4 Handling non-equivalent preprocessing

The Cardinal and MALDIquant preprocessing pipelines may produce different aligned peak sets.

Therefore, perform two separate comparisons.

### Comparison A: Extraction equivalence

Use equivalent numerical spectra and mass axes.

Evaluate the difference caused only by the extraction implementation.

### Comparison B: Complete workflow comparison

Allow the actual MALDIquant and Cardinal preprocessing pipelines to run.

Compare the resulting scientific outputs.

Do not interpret differences from preprocessing as evidence that the analyte extraction implementation is incorrect.

## 25.5 Spatial reference reproducibility

Fix random seeds.

Where possible, use identical input point patterns for the spatial statistical validation.

This isolates the MSI backend migration from downstream stochastic behavior.

---

# 26. Required Test Matrix

The following test matrix must be implemented progressively.

| Component | Required tests |
|---|---|
| Import | Continuous/processed, valid/invalid files, memory modes |
| Data structures | Arrays, shared-axis experiments, metadata |
| Coordinates | Ordering, orientation, missing pixels, multiple runs |
| Ion extraction | Gaussian, sum, mean, max |
| Mass matching | Exact, approximate, boundary, no match |
| Verification | METASPACE matching, confirmed-only behavior |
| FWHM | Default models, measured widths, interpolation |
| Normalization | TIC, RMS, median, IS, median fold-change |
| Filtering | Peak frequencies and threshold boundaries |
| Batch search | Correctness, ordering, duplicate targets |
| Lipid search | Adducts, database metadata, annotations |
| Point patterns | Coordinates, marks, zero-intensity handling |
| MPMs | Density estimates and spatial classifications |
| Shiny | Import, initialization, analysis, session lifecycle |
| Memory | Large file-backed datasets |
| Parallelization | Serial and parallel equivalence |
| Packaging | Installation, checks, dependency validation |

## 26.1 Essential edge cases

Explicitly test:

- Empty mass spectra.
- All-zero spectra.
- Very large intensity values.
- Very small intensity values.
- `NA`, `NaN`, and `Inf`.
- Duplicate peaks.
- Unsorted mass values.
- Nonuniform m/z grids.
- Missing pixels.
- Duplicate coordinates.
- Multiple runs.
- Different pixel spacing.
- Out-of-range analytes.
- Peaks exactly at tolerance boundaries.
- Zero or negative FWHM.
- Extremely narrow search windows.
- Extremely broad search windows.
- Datasets with only one pixel.
- Datasets with only one detected feature.

The implementation must either handle these cases correctly or reject them with informative errors.

---

# 27. Performance and Scalability Requirements

## 27.1 Primary performance objective

The new package must support MSI datasets that are impractical to process through the current MALDIquant-based workflow because of memory limitations.

This must be demonstrated empirically.

## 27.2 Benchmark sizes

Use representative datasets containing approximately:

| Pixel count | Objective |
|---:|---|
| 5,000 | Basic correctness and overhead |
| 50,000 | Moderate dataset |
| 100,000 | Large dataset |
| 200,000 | Large-scale performance |
| 500,000+ | Stress testing |

These are proposed benchmark targets.

They are not guaranteed performance capabilities.

## 27.3 Required measurements

For each benchmark, record:

```text
Dataset size
Number of pixels
Average peaks per spectrum
Total spectral measurements
File size
Import time
Preprocessing time
Single-analyte extraction time
Batch extraction time
MPM computation time
Peak RAM consumption
Disk usage
Operating system
R version
Cardinal version
matter version
CPU model
Available RAM
Storage type
```

## 27.4 Memory-efficiency criteria

Avoid operations that:

- Materialize complete spectral matrices unnecessarily.
- Create global vectors containing all peak positions.
- Repeatedly copy large intensity arrays.
- Read the same dataset once for every analyte in a large batch.
- Duplicate the complete dataset in parallel workers.

A large dataset should remain analyzable without requiring sufficient RAM to hold the entire spectral dataset.

## 27.5 Expected limitations

The new architecture may still encounter substantial computation costs in:

- Spatial kernel density estimation.
- CSR simulations.
- Very large collective projections.
- Plotting high-resolution spatial results.
- Extremely large batch queries.

These costs should be measured separately from spectral import and extraction.

Do not attribute every runtime limitation to Cardinal.

---

# 28. Implementation Roadmap

The development must follow the sequence below.

## Phase 0 — Repository preparation

### Tasks

- Confirm the working branch is `Cardinal-dev`.
- Preserve the legacy implementation.
- Record the starting commit.
- Record dependency versions.
- Establish a reproducible legacy environment.
- Identify representative test datasets.

### Deliverables

- Repository audit.
- Dependency audit.
- Legacy validation fixtures.
- Initial migration documentation.

### Completion criterion

The original implementation can be executed reproducibly and produces fixed reference results.

---

## Phase 1 — Cardinal infrastructure

### Tasks

- Introduce Cardinal as a dependency.
- Implement `readMsi()`.
- Support file-backed import.
- Implement coordinate and metadata access helpers.
- Validate `MSImagingArrays`.
- Validate `MSImagingExperiment`.
- Support run selection.

### Deliverables

- Cardinal import implementation.
- Input validation.
- Unit tests.
- Basic documentation.

### Completion criterion

moleculaR can import and inspect MSI data using Cardinal without constructing MALDIquant objects.

---

## Phase 2 — Single-analyte extraction

### Tasks

- Implement Cardinal-native intensity extraction.
- Support Gaussian weighting.
- Support sum, mean, and max aggregation.
- Preserve FWHM-based search.
- Preserve annotation verification.
- Preserve point-pattern metadata.
- Support both Cardinal representations.
- Validate spatial coordinates.

### Deliverables

- New `searchAnalyte()` implementation.
- Numerical regression tests.
- Small working MPM example.

### Completion criterion

A Cardinal MSI object can be passed directly through:

```text
searchAnalyte()
      |
      v
probMap()
```

without using MALDIquant or creating a global sparse intensity matrix.

**This is the first major project milestone.**

---

## Phase 3 — Batch extraction

### Tasks

- Implement chunked multi-analyte extraction.
- Refactor `batchAnalyteSearch()`.
- Refactor `batchLipidSearch()`.
- Introduce optional caching.
- Benchmark disk access.
- Evaluate parallel processing.

### Deliverables

- Batch extraction engine.
- Updated batch APIs.
- Batch validation tests.
- Initial performance report.

### Completion criterion

Batch extraction produces results equivalent to repeated single-analyte extraction and demonstrates an appropriate memory profile on file-backed datasets.

---

## Phase 4 — Preprocessing and FWHM

### Tasks

- Replace MALDIquant-dependent normalization.
- Replace MALDIquant-dependent peak filtering.
- Adapt `estimateFwhm()`.
- Rewrite `readSingleSpect()`.
- Preserve default FWHM models.
- Validate optional Cardinal peak processing.
- Preserve preprocessing metadata.

### Deliverables

- Cardinal-compatible preprocessing.
- Updated FWHM implementation.
- Numerical regression tests.
- Documentation.

### Completion criterion

The main scientific processing workflow no longer requires MALDIquant.

---

## Phase 5 — Shiny integration

### Tasks

- Replace MALDIquant imports.
- Replace sparse matrix initialization.
- Update mass selection.
- Update preprocessing controls.
- Support file-backed data.
- Implement file lifecycle management.
- Add caching.
- Update example application.

### Deliverables

- Functional Cardinal-based Shiny application.
- Automated tests where practical.
- User documentation.

### Completion criterion

Users can upload supported imzML data and generate MPMs and CPPMs through the Shiny application without MALDIquant.

---

## Phase 6 — Cleanup and documentation

### Tasks

- Remove MALDIquant dependencies.
- Remove unused compatibility code.
- Update package documentation.
- Update examples.
- Update vignettes.
- Update example datasets.
- Update `DESCRIPTION`.
- Regenerate `NAMESPACE`.
- Update `renv.lock`.
- Run package checks.

### Deliverables

- Clean Cardinal-based package.
- Updated installation instructions.
- Updated user guide.
- Migration notes.
- Final validation report.

### Completion criterion

The package installs and passes its tests without MALDIquant.

---

## Phase 7 — Large-scale validation

### Tasks

- Benchmark representative datasets.
- Measure memory usage.
- Compare single and batch extraction.
- Assess disk I/O.
- Evaluate chunk sizes.
- Evaluate parallel processing.
- Test large Shiny workflows.

### Deliverables

- Performance benchmark report.
- Scalability assessment.
- Known limitations.
- Recommended computational settings.

### Completion criterion

The migration demonstrates a meaningful improvement in the ability to analyze large MSI datasets.

---

# 29. Suggested Repository Organization

The exact file organization may be adjusted after inspecting the current package structure.

A possible arrangement is:

```text
moleculaR/
│
├── R/
│   ├── readMsi.R
│   ├── msi-accessors.R
│   ├── msi-validation.R
│   ├── msi-extractIon.R
│   ├── msi-extractIons.R
│   ├── msi-processing.R
│   ├── searchAnalyte.R
│   ├── batchAnalyteSearch.R
│   ├── batchLipidSearch.R
│   ├── estimateFwhm.R
│   ├── useDefaultFwhm.R
│   ├── gaussWeight.R
│   ├── probMap.R
│   └── ...
│
├── tests/
│   ├── testthat/
│   │   ├── test-cardinal-import.R
│   │   ├── test-cardinal-arrays.R
│   │   ├── test-cardinal-experiment.R
│   │   ├── test-analyte-extraction.R
│   │   ├── test-batch-extraction.R
│   │   ├── test-normalization.R
│   │   ├── test-fwhm.R
│   │   ├── test-spatial-integration.R
│   │   └── test-file-backed.R
│   │
│   └── fixtures/
│       └── legacy/
│
├── benchmarks/
│   ├── benchmark-import.R
│   ├── benchmark-extraction.R
│   ├── benchmark-batch.R
│   └── benchmark-memory.R
│
├── docs/
│   └── Cardinal-migration/
│       ├── dependency-audit.md
│       ├── architecture.md
│       ├── scientific-decisions.md
│       ├── numerical-validation.md
│       ├── performance.md
│       └── migration-report.md
│
├── vignettes/
│   └── ...
│
├── DESCRIPTION
├── NAMESPACE
└── CARDINAL_MIGRATION_AGENT.md
```

Avoid creating unnecessary abstraction layers.

The structure should remain understandable and maintainable by researchers familiar with R and MSI analysis.

---

# 30. Detailed Risk Register

## Risk 1 — Different peak-alignment results

**Cause:** Cardinal and MALDIquant use different peak alignment and binning algorithms.

**Potential consequence:** Different feature positions and intensities.

**Mitigation:** Avoid mandatory global alignment. Validate preprocessing independently from extraction.

## Risk 2 — Different mass-window selection

**Cause:** Legacy nearest-index matching versus strict numerical tolerance.

**Potential consequence:** Different analyte intensities.

**Mitigation:** Unit tests with controlled boundary cases; document the chosen rule.

## Risk 3 — Different mean aggregation

**Cause:** Legacy global-feature denominator versus per-spectrum observed peaks.

**Potential consequence:** Different normalized mean intensities.

**Mitigation:** Define and test the intended denominator explicitly.

## Risk 4 — Pixel-coordinate mismatch

**Cause:** Different indexing or matrix orientation.

**Potential consequence:** Incorrect molecular maps.

**Mitigation:** Asymmetric synthetic spatial fixtures and coordinate-level comparisons.

## Risk 5 — Multiple experimental runs

**Cause:** Identical coordinates in different runs.

**Potential consequence:** Incorrect aggregation of separate samples.

**Mitigation:** Explicit run selection and validation.

## Risk 6 — Loss of file-backed behavior

**Cause:** Unnecessary conversion to ordinary R matrices or lists.

**Potential consequence:** High RAM consumption.

**Mitigation:** Inspect data-access paths and benchmark memory.

## Risk 7 — Repeated disk reads

**Cause:** Naive single-analyte extraction repeated thousands of times.

**Potential consequence:** Slow batch analysis.

**Mitigation:** Chunked multi-analyte extraction and caching.

## Risk 8 — Changed FWHM estimates

**Cause:** Different peak-detection algorithms or interpolation.

**Potential consequence:** Different analyte detection.

**Mitigation:** Preserve numerical width calculations and validate fitted curves.

## Risk 9 — Lost example-data FWHM reference

**Cause:** Removing MALDIquant-based bundled data.

**Potential consequence:** Changed default FTICR FWHM model.

**Mitigation:** Extract and preserve the required numerical reference.

## Risk 10 — Deferred preprocessing inconsistencies

**Cause:** Processing steps may remain queued.

**Potential consequence:** Analysis of unintended spectral intensities.

**Mitigation:** Explicit processing-state handling and metadata recording.

## Risk 11 — Shiny temporary-file deletion

**Cause:** File-backed objects reference temporary uploaded files.

**Potential consequence:** Runtime errors after initialization.

**Mitigation:** Managed session storage.

## Risk 12 — Dependency incompatibilities

**Cause:** Mixing incompatible R and Bioconductor versions.

**Potential consequence:** Installation or runtime failure.

**Mitigation:** Pin and validate a compatible environment.

## Risk 13 — Sparse matrix orientation

**Cause:** Cardinal's feature-by-pixel orientation differs from moleculaR's legacy pixel-by-feature representation.

**Potential consequence:** Incorrect indexing or intensity assignment.

**Mitigation:** Explicit orientation tests for every shared-axis operation.

## Risk 14 — Memory use in spatial statistics

**Cause:** Large spatstat point patterns, raster images, or simulations.

**Potential consequence:** Continued performance limitations after spectral migration.

**Mitigation:** Benchmark spectral and spatial processing separately.

---

# 31. Acceptance Criteria

The migration is complete only when all mandatory criteria are satisfied.

## 31.1 Architectural criteria

- [ ] Cardinal is the primary MSI infrastructure.
- [ ] `MSImagingArrays` is supported.
- [ ] `MSImagingExperiment` is supported.
- [ ] File-backed import works.
- [ ] A global sparse matrix is not required for single-analyte extraction.
- [ ] MALDIquant is no longer required for normal operation.
- [ ] MALDIquantForeign is no longer required for normal operation.
- [ ] Cardinal objects remain usable with Cardinal's own API.

## 31.2 Scientific criteria

- [ ] Gaussian weighting is preserved.
- [ ] FWHM-dependent mass matching is preserved.
- [ ] Alternative aggregation methods are supported.
- [ ] Annotation verification is preserved.
- [ ] Point-pattern structure is preserved.
- [ ] MPM generation remains functional.
- [ ] CPPM generation remains functional.
- [ ] Lipid metadata are preserved.
- [ ] Scientific differences from legacy behavior are documented.

## 31.3 Technical criteria

- [ ] Import tests pass.
- [ ] Extraction tests pass.
- [ ] Batch extraction tests pass.
- [ ] FWHM tests pass.
- [ ] Normalization tests pass.
- [ ] Spatial integration tests pass.
- [ ] Shiny workflows pass.
- [ ] R package checks pass.
- [ ] Supported platforms are tested.
- [ ] Reproducible dependency configuration is provided.

## 31.4 Performance criteria

- [ ] Large file-backed datasets can be imported without full spectral materialization.
- [ ] Single-analyte extraction does not require loading the complete dataset into RAM.
- [ ] Batch extraction avoids unnecessary repeated full-dataset reads.
- [ ] Memory usage is measured and documented.
- [ ] Large-scale benchmarks are completed.
- [ ] Known scalability limitations are documented.

---

# 32. Definition of Done

The project is considered complete when a user can execute a workflow equivalent to:

```r
library(moleculaR)
library(Cardinal)

msData <- readMsi(
    path = "dataset.imzML",
    memory = FALSE
)

fwhmObj <- useDefaultFwhm(
    modality = "FTICR"
)

mzTarget <- 788.5447

analyte <- searchAnalyte(
    m = mzTarget,
    fwhm = getFwhm(fwhmObj, mzTarget),
    spData = msData
)

mpm <- probMap(
    analyte,
    seed = 1
)

plot(mpm, what = "MPM")
```

without:

- Importing MALDIquant.
- Converting the dataset into MALDIquant spectra.
- Running `MALDIquant::binPeaks()`.
- Constructing a global sparse intensity matrix.

The package must also support large batch searches and the existing molecular probability mapping functionality.

The numerical behavior must be validated against the legacy implementation.

---

# 33. Instructions for the First Development Session

The agent should begin with the following tasks.

## Step 1

Inspect the complete repository.

Identify the current MALDIquant dependencies.

Create the dependency audit.

## Step 2

Inspect the Cardinal documentation.

Confirm the installed Cardinal, matter, and R versions.

Verify the exact APIs for:

```r
Cardinal::readMSIData()
Cardinal::MSImagingArrays()
Cardinal::MSImagingExperiment()
Cardinal::mz()
Cardinal::intensity()
Cardinal::coord()
Cardinal::spectrapply()
Cardinal::process()
```

## Step 3

Create a small synthetic `MSImagingArrays` object.

Confirm that its spectra and coordinates can be retrieved.

## Step 4

Implement a prototype of direct Gaussian-weighted analyte extraction.

Do not modify `probMap()`.

## Step 5

Convert the extracted intensity vector to the existing `analytePointPattern` class.

## Step 6

Generate an MPM from the Cardinal-derived point pattern.

## Step 7

Compare the extracted intensities and MPM with the legacy implementation on a small controlled dataset.

## Step 8

Document the results.

Only after completing these steps should the agent proceed to batch extraction and the remaining migration tasks.

---

# 34. Final Architectural Guidance

The most important principle of the migration is:

**Cardinal should manage the MSI data. moleculaR should manage the molecular probability mapping methodology.**

Avoid reproducing the MALDIquant architecture inside Cardinal.

The central improvement comes from eliminating unnecessary global spectral conversion and allowing moleculaR to operate directly on file-backed MSI datasets.

The migration should therefore prioritize:

1. Correct Cardinal data access.
2. Scientifically validated analyte extraction.
3. Efficient batch processing.
4. Preservation of spatial analysis behavior.
5. Long-term maintainability.

The first and most important milestone is a fully functional Cardinal-native implementation of `searchAnalyte()`.

Once this is achieved and validated, the remaining migration can proceed incrementally.

**The development agent must prioritize scientific correctness and memory-efficient architecture over superficial API compatibility or premature optimization.**

The final result should be a modernized moleculaR package capable of analyzing large MSI datasets while retaining the scientific framework established by the original implementation.

---

# 35. Required Final Deliverables

At the end of the migration, provide:

1. A functional Cardinal-based moleculaR package.
2. Updated source code and dependencies.
3. Updated documentation.
4. Updated vignettes.
5. Updated example workflows.
6. Updated Shiny applications.
7. Automated unit and integration tests.
8. Legacy numerical comparison results.
9. A performance benchmarking report.
10. A migration report describing all major changes.
11. A list of known limitations.
12. A reproducible installation environment.

The migration report must explicitly state:

- Which functionality is scientifically equivalent to the legacy implementation.
- Which functionality produces different numerical results.
- Why those differences occur.
- Whether the differences are expected.
- Whether they affect molecular probability maps.
- How the package handles very large datasets.

Do not mark the project complete with unresolved scientific discrepancies.

---

**End of development specification.**