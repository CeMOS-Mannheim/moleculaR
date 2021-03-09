## Classes
##
## Under the hood, moleculaR relies on MALDIquant, MALDIquantForeign,
## Matrix and spatstat packages.
##
## As of the current version, moleculaR does not declare any new classes,
## instead it relies on 'MALDIquant::MassPeaks' class to internally represent
## centroided mass spectra. moleculaR does not support continuous MSI data,
## instead, it assumes centroided MSI data as an input.
##
## It also uses class definitions of 'Matrix' and 'spatstat' packages. For
## more info please check the documentation of the respective packages.
##
##
