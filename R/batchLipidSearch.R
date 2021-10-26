#' Batch Analyte detection in an MSI dataset
#'
#' This applies `moleculaR::searchAnalyte` on a given MSI dataset against the entire SwissLipids database
#' taking into account the different ionization statuses (i.e. adducts).
#'
#' @param spData an S3 object of type `sparseIntensityMatrix` holding the sparse MSI data.
#' @param fwhmObj an S3 object of type `fwhm` with the estimated fwhm data.
#' @param sldb the SwissLipid database loaded as a data.frame.
#' @param adduct a character vector specifying the ionisation status of interest, c("M-H", "M+H", "M+Na", "M+K").
#' @param numCores an integer, the number of cores to be used for the search, passed to `parallel::mclapply`. Not
#' supported on Windows machines.
#' @param wMethod wighting method; c("Gaussian", "sum", "max", "mean").
#' @param verifiedMasses an optional numeric vector of m/z values that are (externally) verified to
#' be real molecular entities with a certain confidence, ex. 'mz' column of a METASPACE
#' annotation result.
#' @param confirmedOnly if `TRUE`, returns detections only if confirmed by `verifiedMasses`.
#' @param metaData optional named list with additional identifiers for the analyte
#' under study, ex. list(lipidID = "..", sumformular = "..", ..). This will be passed to the `metaData`
#' slot of the resulting `analytePointPattern` object.
#'
#' @return An analyte point patter of type `ppp` and `analytePointPattern` containing all lipid hits identified in the
#' MSI dataset `spData` for the specified `adduct` formations taking into account the `fwhm` information.
#'
#' @export
#'

batchLipidSearch <- function(spData, fwhmObj, sldb, adduct = c("M-H", "M+H", "M+Na", "M+K"),
                             numCores = 1L, wMethod = "Gaussian", verifiedMasses = NA,
                             confirmedOnly = FALSE, verbose = TRUE) {

      #>>> to do: confirmedOnly flag set manually here. The check has to be moved to search Analyte.


      #// sort verifiedMasses
      verifiedMasses <- sort(verifiedMasses)

      #// create sp window
      spwin <- spatstat::as.polygonal(spatstat::owin(mask = spData$coordinates))


      #// filter sldb to include only verified masses - speed up computations
      if(confirmedOnly){

            if(identical(verifiedMasses, NA)){
                  stop("'verifiedMasses' has to be provided when 'confirmedOnly'  is TRUE.\n")
            }

            Hidx <- MALDIquant::match.closest(x = sldb$`Exact m/z of [M+H]+`,
                                             table = sort(verifiedMasses),
                                             tolerance = getFwhm(fwhmObj, sldb$`Exact m/z of [M+H]+`))

            Naidx <- MALDIquant::match.closest(x = sldb$`Exact m/z of [M+Na]+`,
                                              table = sort(verifiedMasses),
                                              tolerance = getFwhm(fwhmObj, sldb$`Exact m/z of [M+Na]+`))

            kidx <- MALDIquant::match.closest(x = sldb$`Exact m/z of [M+K]+`,
                                             table = sort(verifiedMasses),
                                             tolerance = getFwhm(fwhmObj, sldb$`Exact m/z of [M+K]+`))

            idxToKeep <- c(which(!is.na(Hidx)), which(!is.na(Naidx)), which(!is.na(kidx)))

            sldb <- sldb[idxToKeep, ]

      }


      #// initialize the swisslipids database
      searchList <- initLipidSearch(swissdb = sldb)


      # ofInterest = c("FA(x:x)", "LPA(x:x)", "LPC(x:x)", "LPE(x:x)", "LPG(x:x)",
      #                "LPI(x:x)", "LPS(x:x)", "PA(x:x)", "PC(x:x)", "PE(x:x)",
      #                "PG(x:x)", "PI(x:x)", "PS(x:x)", "DG(x:x)", "TG(x:x)",
      #                "PGP(x:x)", "PIP2(x:x)", "PIP(x:x)", "PIP3(x:x)")

      if(verbose){
            pb <- utils::txtProgressBar(min = 1, max = length(names(searchList$hitsList)), style = 3, width = 10)
            count <- 1
      }

      for(lipidClass in names(searchList$hitsList)){

            if(verbose)
                  utils::setTxtProgressBar(pb, count); count = count +1

            #if(!(lipidClass %in% ofInterest)) {next}

            searchList$hitsList[[lipidClass]] <- parallel::mclapply(X = seq(1, nrow(searchList$lipidList[[lipidClass]])),
                                                                    mc.cores = numCores, FUN = function(i) {


                          sppCotainer <- list(emptyspp = spatstat::ppp(x = integer(0), y = integer(0)))

                          #// de-protonated ----
                          if("M-H" %in% adduct) {

                                lipTmp    <- searchList$lipidList[[lipidClass]]$`Exact m/z of [M-H]-`[i]

                                if(!is.na(lipTmp)) { # for example there is no protonated version of this lipid species

                                      # metaData list
                                      mt  <- list(mode = "negative",
                                                  adduct = "M-H",
                                                  lipidID = searchList$lipidList[[lipidClass]]$`Lipid ID`[i],
                                                  sumformula = searchList$lipidList[[lipidClass]]$`Formula (pH7.3)`[i],
                                                  abbrev = searchList$lipidList[[lipidClass]]$`Abbreviation*`[i],
                                                  numDoubleBonds = searchList$lipidList[[lipidClass]]$numDoubleBond[i])

                                      hitsDeprot <- searchAnalyte(m = lipTmp,
                                                                fwhm = getFwhm(fwhmObj, lipTmp),
                                                                spData = spData,
                                                                wMethod = wMethod,
                                                                verifiedMasses = verifiedMasses,
                                                                confirmedOnly = confirmedOnly,
                                                                metaData = mt)

                                      sppCotainer <- c(sppCotainer, list(hitsDeprot))

                                      rm(mt, lipTmp)


                                }

                          }


                          #// protonated ----
                          if("M+H" %in% adduct){
                                lipTmp        = searchList$lipidList[[lipidClass]]$`Exact m/z of [M+H]+`[i]

                                if(!is.na(lipTmp)) { # for example there is no protonated version of this lipid species


                                      # metaData list
                                      mt  <- list(mode = "positive",
                                                  adduct = "M+H",
                                                  lipidID = searchList$lipidList[[lipidClass]]$`Lipid ID`[i],
                                                  sumformula = searchList$lipidList[[lipidClass]]$`Formula (pH7.3)`[i],
                                                  abbrev = searchList$lipidList[[lipidClass]]$`Abbreviation*`[i],
                                                  numDoubleBonds = searchList$lipidList[[lipidClass]]$numDoubleBond[i])

                                      hitsProt <- searchAnalyte(m = lipTmp,
                                                              fwhm = getFwhm(fwhmObj, lipTmp),
                                                              spData = spData,
                                                              wMethod = wMethod,
                                                              verifiedMasses = verifiedMasses,
                                                              confirmedOnly = confirmedOnly,
                                                              metaData = mt)

                                      sppCotainer <- c(sppCotainer, list(hitsProt))


                                      rm(mt, lipTmp)



                                }
                          }



                          #// Na+ adduct ----
                          if("M+Na" %in% adduct){
                                lipTmp               = searchList$lipidList[[lipidClass]]$`Exact m/z of [M+Na]+`[i]

                                if(!is.na(lipTmp)) { # for example there is no Na-adduct version



                                      # metaData list
                                      mt  <- list(mode = "positive",
                                                  adduct = "M+Na",
                                                  lipidID = searchList$lipidList[[lipidClass]]$`Lipid ID`[i],
                                                  sumformula = searchList$lipidList[[lipidClass]]$`Formula (pH7.3)`[i],
                                                  abbrev = searchList$lipidList[[lipidClass]]$`Abbreviation*`[i],
                                                  numDoubleBonds = searchList$lipidList[[lipidClass]]$numDoubleBond[i])

                                      hitsSod <- searchAnalyte(m = lipTmp,
                                                             fwhm = getFwhm(fwhmObj, lipTmp),
                                                             spData = spData,
                                                             wMethod = wMethod,
                                                             verifiedMasses = verifiedMasses,
                                                             confirmedOnly = confirmedOnly,
                                                             metaData = mt)

                                      sppCotainer <- c(sppCotainer, list(hitsSod))


                                      rm(mt, lipTmp)



                                }

                          }


                          #// K+ adduct ----
                          if("M+K" %in% adduct | "M+k" %in% adduct){
                                lipTmp        = searchList$lipidList[[lipidClass]]$`Exact m/z of [M+K]+`[i]

                                if(!is.na(lipTmp)) { # for example there is no Na-adduct version

                                      # metaData list
                                      mt  <- list(mode = "positive",
                                                  adduct = "M+K",
                                                  lipidID = searchList$lipidList[[lipidClass]]$`Lipid ID`[i],
                                                  sumformula = searchList$lipidList[[lipidClass]]$`Formula (pH7.3)`[i],
                                                  abbrev = searchList$lipidList[[lipidClass]]$`Abbreviation*`[i],
                                                  numDoubleBonds = searchList$lipidList[[lipidClass]]$numDoubleBond[i])

                                      hitsPotas <- searchAnalyte(m = lipTmp,
                                                               fwhm = getFwhm(fwhmObj, lipTmp),
                                                               spData = spData,
                                                               wMethod = wMethod,
                                                               verifiedMasses = verifiedMasses,
                                                               confirmedOnly = confirmedOnly,
                                                               metaData = mt)

                                      sppCotainer <- c(sppCotainer, list(hitsPotas))


                                      rm(mt, lipTmp)



                                }
                          }

                          if(length(sppCotainer) > 1) {
                             sppCotainer <- superImposeAnalytes(sppCotainer, spWin = spwin)
                          }

                          return(sppCotainer)
                        })





            #// merge
            searchList$hitsList[[lipidClass]] <- superImposeAnalytes(searchList$hitsList[[lipidClass]],
                                                                      spWin = spwin)

            if(confirmedOnly){ # manually set this flag
                  searchList$hitsList[[lipidClass]]$metaData$mzConfirmed <- TRUE
            }
      }

      if(verbose)
            close(pb)

      searchList <- .reduceSearchList(searchList)

      return(searchList)


}


.reduceSearchList <- function(searchList){ # remove classes which did not generate hits

      tokeep <- !(sapply(searchList$hitsList, function(x) spatstat::is.empty.ppp(x)))

      searchList <- lipidSearchList(lipidList = searchList$lipidList,
                                    hitsList = searchList$hitsList[tokeep],
                                    allClasses = names(searchList$hitsList)[tokeep])

      return(searchList)

}
