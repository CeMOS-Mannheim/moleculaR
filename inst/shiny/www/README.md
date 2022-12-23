

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

The first step involves uploading data. To do this click on `imzML & ibd` File and navigate to the directory containing the files you want to upload. Then select both imzML- and ibd-Files and upload them. Via `Spectrum .tsv File` a continuous spectrum in `tsv` format could be uploaded. This spectrum could either be a single random pixel or a mean spectrum of your imaging dataset. Once this is done, click `Load and Initialize` to read the data and estimate the FWHM model. For a given m/z value of interest, you can generate the corresponding MPM by providing that m/z value in the corresponding input box. To generate a collective projection probabilstic map (CPPM) of a custom list of m/z values, please paste these values into the corresponding text box, comma separated.


### Citing _moleculaR_

Abu Sammour, Denis, et al. "Spatial Probabilistic Mapping of Metabolite Ensembles in Mass Spectrometry Imaging." bioRxiv (2021). https://doi.org/10.1101/2021.10.27.466114 


### Contact

You are welcome to:

* submit suggestions and bug-reports at: <https://github.com/CeMOS-Mannheim/moleculaR/issues>
* send a pull request on: <https://github.com/CeMOS-Mannheim/moleculaR/>
* compose an e-mail to: <d.abu-sammour@hs-mannheim.de>



