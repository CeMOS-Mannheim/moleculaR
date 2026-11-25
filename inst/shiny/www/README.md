

<p align="center"><img src="moleculaR-icon1.svg" width="400"></p>


<br />
<br />

The _moleculaR_ R package provides a computational framework that introduces probabilistic mapping
and point-for-point statistical testing of metabolites in tissue via Mass spectrometry imaging.
It enables collective projections of metabolites and consequently spatially-resolved investigation
of ion milieus, lipid pathways or user-defined biomolecular ensembles within the same image.

For more information about this package and its applications please refer to the [associated preprint](https://doi.org/10.1101/2021.10.27.466114). 



### _moleculaR_ Shiny-App

_moleculaR_ provides an R [Shiny](https://www.rdocumentation.org/packages/Shiny/versions/1.7.1) web app `package-app` with intuitive web-based GUI. It contains an `example` section which comes pre-loaded with an examplary reduced MALDI MSI dataset (see citation below) and a `main` section which lets the user upload her own centroided [imzML](https://ms-imaging.org/imzml/) data and apply spatial probabilistic mapping through Molecular probabilistic Maps (MPMs) and Collective Projection probabilistic Maps (CPPMs). 

Both apps use native, file-backed Cardinal data. Upload one matching **known-centroided** imzML/ibd pair, select an experimental run and optional coordinate ROI, then click **Load and initialize**. Profile or unknown-status MSI must be explicitly processed in Cardinal before upload. All acquired pixels, including zero-signal pixels and tissue holes, define the selected domain. Supply pixel spacing if the grid cannot be inferred.

Under **FWHM source**, select **Approximate instrument model** and FTICR, TOF-reflector, or TIMS-TOF. Default resolving powers are 123000, 10000, and 40000 at m/z 249; custom calibration is optional. For measured widths, select **Separate profile spectrum** and upload a two-column profile TSV (m/z, intensity), explicitly choosing whether it has a header. Cardinal picks peaks only in that separate profile, using the displayed picker and SNR. All valid selected widths are fitted without random subsampling. Reduced centroided imzML, including the bundled example, is unsuitable for peak picking or width estimation. Failed estimation reports an error.

Normalization defaults to **None**. Optional Cardinal TIC, RMS, median or internal-standard normalization creates a private intensity file and preserves the source. Reinitialization starts from the original copied input. TIC on reduced data sums only the supplied features. Changing data, ROI, normalization, annotations or FWHM settings requires initialization again.

Enter an exact theoretical m/z for an MPM, or a comma/space-separated mass list for a CPPM. Masses are never snapped to an acquired feature. Choose Gaussian, sum, max or mean aggregation; optional annotation confirmation; bandwidth; and a random seed (default 1). Collective maps can standardize each analyte using the existing z-score/range transform. Constant or insufficient signals require disabling standardization or revising the query. Query changes clear results; changing the display only redraws the map.

The Lipids tab supports class, (lyso)GPL ion and saturation groups. Load the bundled SwissLipids species metadata or an uploaded SwissLipids TSV; the example uses a smaller annotation-matched subset. Lipid identifiers/adducts remain separate at repeated masses. Annotation CSVs need an `mz` column and must already be filtered to your intended FDR; METASPACE exports with two metadata rows are supported.

Download current results as a portable analysis RDS, an intensity/metadata CSV or a probability-map PNG. The RDS contains numeric spatial results, FWHM and provenance, without the live MSI object. Each tab owns its copied uploads, cache and derived files; session closure removes those private files. Computation is synchronous within each R worker, with staged progress. Configure `options(shiny.maxRequestSize = ...)` before launch to choose an upload limit; the app supplies a 48 GiB default only when that option is unset.


### Citing _moleculaR_

Abu Sammour, Denis, et al. "Spatial Probabilistic Mapping of Metabolite Ensembles in Mass Spectrometry Imaging." bioRxiv (2021). https://doi.org/10.1101/2021.10.27.466114 


### Contact

You are welcome to:

* submit suggestions and bug-reports at: <https://github.com/CeMOS-Mannheim/moleculaR/issues>
* send a pull request on: <https://github.com/CeMOS-Mannheim/moleculaR/>
* compose an e-mail to: <d.abu-sammour@hs-mannheim.de>



