This example is a 32-by-32 coordinate crop (1,024 acquired pixels) of the preserved
reduced centroided glioblastoma example. It is unsuitable for peak picking or
measured FWHM estimation. The apps use approximate instrument calibration by
default; a separate profile TSV is required for measured widths.

`example.imzML` and `example.ibd` store native processed centroid spectra with
float64 intensities. Keep the pair together. `annotations.csv` contains the
original verified masses/FDR metadata. `lipids.rds` contains 546 species with
positive-adduct theoretical masses matching those annotations within three
default-FTICR sigmas. `swisslipids-species.rds` contains all 10,319 species from
the original bundled database for Main's default database. Both are ordinary
data frames with the same ten mass/identifier columns. These metadata are a
preserved reference, not a newly downloaded database release.

`provenance.rds` records the original data checksum, crop bounds, source commit
and processing assumptions. No normalization, picking or width estimation was
performed when exporting the crop. The development-only generator is
`docs/Cardinal-migration/generate-shiny-example.R`; app runtime never loads the
legacy spectral bundle. The original dataset is described in the package's
glioblastoma example documentation and Abu Sammour et al., Nature Communications
(2023), https://doi.org/10.1038/s41467-023-37394-z.
