# **moleculaR**

### Spatial Probabilistic Mapping of Metabolite Ensembles in Mass Spectrometry Imaging

<br />
<br />

The moleculaR R package provides a computational framework that introduces probabilistic mapping
and point-for-point statistical testing of metabolites in tissue via Mass spectrometry imaging.
It enables collective molecular projections and consequently spatially-resolved investigation
of ion milieus, lipid pathways or user-defined biomolecular ensembles within the same image.

MoleculaR comes pre-loaded with the [SwissLipids database](https://www.swisslipids.org) and with is capable of importing metabolite annotation results from the [METASPACE platform](https://metaspace2020.eu/) to compute FDR-verified molecular probability maps (MPMs) and collective projection probability maps (CPPMs). moleculaR could also be deployed and hosted on a centralized server and is equipped with a web-based GUI based on [shiny](https://www.rdocumentation.org/packages/shiny/versions/1.7.1). 

For more information about this package and its applications please refer to the associated paper (doi). 

<p align="right"><img src="extras/package.jpg" width="680"></p>

### Installation

The [devtools](https://cran.r-project.org/web/packages/devtools/index.html) package could be used to install the development version of moleculaR:

```r
install.packages("devtools")

library("devtools")
install_github("CeMOS-Mannheim/moleculaR")
```
Note that `moleculaR` was created with `renv`to help manage R package dependencies and computational reproducibility. for more info, check [the renv guide page](https://rstudio.github.io/renv/articles/renv.html). 



### Vignettes



### running shiny example-app



### running shiny package-app


### deplyoing on premises

`moleculaR` comes pre-loaded with two R Shiny-apps which could also be hosted on a local server by downloading and installing [Shiny Server](https://www.rstudio.com/products/shiny/shiny-server/) with which a user could run these apps from a local (or possibly remote) network from their browser. For more info please refer to the [Shiny Server download page](https://www.rstudio.com/products/shiny/download-server/) and [installation instruction and server management](https://docs.rstudio.com/shiny-server/#installation). 


### Citing moleculaR



### Contact

You are welcome to:

* submit suggestions and bug-reports at: <https://github.com/CeMOS-Mannheim/moleculaR/issues>
* send a pull request on: <https://github.com/CeMOS-Mannheim/moleculaR/>
* compose an e-mail to: <d.abu-sammour@hs-mannheim.de>

### License

See [license document](LICENSE.md).


