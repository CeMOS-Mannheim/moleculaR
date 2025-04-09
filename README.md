

<p align="center"><img src="extras/moleculaR-icon1.svg" width="400"></p>

### Spatial Probabilistic Mapping of Metabolite Ensembles in Mass Spectrometry Imaging

<br />
<br />

The _moleculaR_ R package provides a computational framework that introduces probabilistic mapping
and point-for-point statistical testing of metabolites in tissue via Mass spectrometry imaging.
It enables collective projections of metabolites and consequently spatially-resolved investigation
of ion milieus, lipid pathways or user-defined biomolecular ensembles within the same image.

_moleculaR_ comes pre-loaded with the [SwissLipids database](https://www.swisslipids.org) and with is capable of importing metabolite annotation results from the [METASPACE platform](https://metaspace2020.eu/) to compute FDR-verified _moleculaR_ probabilistic maps (MPMs) and collective projection probabilistic maps (CPPMs). _moleculaR_ could also be deployed and hosted on a centralized server and is equipped with a web-based GUI based on [Shiny](https://www.rdocumentation.org/packages/Shiny/versions/1.7.1). 

For more information about this package and its applications please refer to the [associated published article](https://doi.org/10.1038/s41467-023-37394-z). 

<p align="center"><img src="extras/package.jpg" width="400"></p>

### Installation

The [devtools](https://cran.r-project.org/web/packages/devtools/index.html) package could be used to install the development version of _moleculaR_. By using `build_vignettes=TRUE`, the vignettes provided with the package will be automatically build, but keep in mind that this will prolong the package installation processes. 

```r
install.packages("devtools")
devtools::install_github("CeMOS-Mannheim/moleculaR", build_vignettes=TRUE)
```
Note that _moleculaR_ was created with `renv`to help manage R package dependencies and computational reproducibility. for more info, check [the renv guide page](https://rstudio.github.io/renv/articles/renv.html). 

### Troubleshooting
#### _XML_ Installation Issues on Windows

The package _MALDIquantForeign_ is a suggested dependency of _moleculaR_, which in turn 
depends on the _XML_ package. On Windows systems, _XML_ may fail to compile from source. 
To avoid this issue, it is recommended to install _XML_ as a binary package:

```r
install.packages("XML", type = "binary")
```

Alternatively, when using _renv_, use the following:

```r
renv::install("XML", type = "binary")
```

Due to this compilation issue, _MALDIquantForeign_ has been moved from the `Imports` to the 
`Suggests` section of the package dependencies.

#### _GDAL_ Dependency on Linux

On Linux-based systems, the R package _terra_ (a dependency of _moleculaR_) requires _GDAL_-related system libraries. 
These may not be pre-installed and could cause installation issues. On Debian-based systems (e.g., Ubuntu), 
you can install the required libraries via:

```r
sudo apt install gdal-bin libgdal-dev
```

Make sure these system dependencies are in place before attempting to install _moleculaR_.

### Data availability

Example MSI data could be downloaded in [imzML](https://ms-imaging.org/imzml/) fomat via [this link](https://metaspace2020.eu/project/abusammour-2021) (not yet public). 

### Vignettes

Vignettes are provided with the package to illustrate basic functionality. You can see all installed vignettes by calling `browseVignettes("moleculaR")`. To read a speficic vignette use, for example, `vignette("moleculaR-walkthrough")` and to view its code use `edit(vignette("moleculaR-walkthrough"))`. 

### _moleculaR_ Shiny-Apps

_moleculaR_ provides an R [Shiny](https://www.rdocumentation.org/packages/Shiny/versions/1.7.1) web app `package-app` with intuitive web-based GUI. It contains an `example` section which comes pre-loaded with an examplary reduced MALDI MSI dataset (see citation below) and a `main` section which lets the user upload her own centroided [imzML](https://ms-imaging.org/imzml/) data and apply spatial probabilistic mapping through Molecular probabilistic Maps (MPMs) and Collective Projection probabilistic Maps (CPPMs). 

### Deplyoing on a Local Server

_moleculaR_ comes pre-loaded with two R Shiny-apps which could also be hosted on a local server by downloading and installing [Shiny Server](https://www.rstudio.com/products/Shiny/Shiny-server/) with which a user could run these apps from a local (or possibly remote) network from their browser. For more info please refer to the [Shiny Server download page](https://www.rstudio.com/products/Shiny/download-server/) and [installation instruction and server management](https://docs.rstudio.com/Shiny-server/#installation). 


### Citing _moleculaR_

Please cite the associated published article [Abu Sammour et al., 2023, Nature Communications](https://doi.org/10.1038/s41467-023-37394-z).


### Contact

You are welcome to:

* submit suggestions and bug-reports at: <https://github.com/CeMOS-Mannheim/moleculaR/issues>
* send a pull request on: <https://github.com/CeMOS-Mannheim/moleculaR/>
* compose an e-mail to: <d.abu-sammour@hs-mannheim.de>

### License

See [license document](LICENSE.md).


