#' Batch Analyte detection in an MSI dataset
#'
#' This applies `moleculaR::searchAnalyte` on a given MSI dataset against the entire SwissLipids database
#' taking into account the different ionization statuses (i.e. adducts).
#'
#' @param spData an S3 object of type `sparseIntensityMatrix` holding the sparse MSI data. This can also be a
#' named list `list(positive=..., negative=...)` to allow for lipid search in positive and negative ionisation
#' mode datasets simultaneously.
#' @param fwhmObj an S3 object of type `fwhm` with the estimated fwhm data.
#' @param spwin optional, an object of type `owin`. If not given the function tries to generate the
#' spatial window out of the coordinates of all points of the dataset stored in `spData` (default behavior).
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
#' @param verbose ignored for now.
#'
#' @return An analyte point patter of type `ppp` and `analytePointPattern` containing all lipid hits identified in the
#' MSI dataset `spData` for the specified `adduct` formations taking into account the `fwhm` information.
#'
#' @export
#'

batchLipidSearch <- function(spData, fwhmObj, spwin = NA, sldb, adduct = c("M-H", "M+H", "M+Na", "M+K"),
                             numCores = 1L, wMethod = "Gaussian", verifiedMasses = NA,
                             confirmedOnly = FALSE, verbose = TRUE) {


        # check for sldb
        if(!("numDoubleBond" %in% colnames(sldb))){
              stop("'numDoubleBond' column is not detecetd in 'sldb'. Please load 'sldb' via 'moleculaR::loadSwissDB' function. \n")
        }
        if(!("numDoubleBond" %in% colnames(sldb))){
            stop("'lipidGroup' column is not detecetd in 'sldb'. Please load 'sldb' via 'moleculaR::loadSwissDB' function. \n")
        }

        # check OS type
        if(.Platform$OS.type == "windows" & numCores > 1){
      	warning("Only single-core operation is supported on windows. \n")
      	numCores <- 1L
        }

      #// sort verifiedMasses
        if(!identical(verifiedMasses, NA)){
              verifiedMasses <- sort(verifiedMasses)

        }

      #// create sp window
      if(identical(spwin, NA)){
            spwin <- spatstat.geom::as.polygonal(spatstat.geom::owin(mask = spData$coordinates))
      }


      #// capitalize "M+k" - to handle input error
      if("M+k" %in% adduct){
            adduct[adduct == "M+k"] = toupper(adduct[adduct == "M+k"])
      }


      #// filter sldb to include only verified masses - speed up computations
      if(confirmedOnly){

            if(identical(verifiedMasses, NA)){
                  stop("'verifiedMasses' has to be provided when 'confirmedOnly'  is TRUE.\n")
            }


            # reduce sldb to only contain verified masses
            sldb <- .filtersldb(sldb, verifiedMasses, adduct, fwhmObj)

      }





      # if(verbose){
      #       pb <- utils::txtProgressBar(min = 1, max = nrow(sldb), style = 3, width = 20)
      #       count <- 1
      # }




      hitsList <- parallel::mclapply(X = seq(1, nrow(sldb)),
                                     mc.cores = numCores, FUN = function(i) {

                    # if(verbose)
                    #       utils::setTxtProgressBar(pb, count); count = count +1




                    #sppCotainer <- list(emptyspp = spatstat.geom::ppp(x = integer(0), y = integer(0)))
                    #sppCotainer <- setNames(object = vector("list", length(adduct)), nm = adduct)

                    # initialize empty spp objects for each adduct
                    sppCotainer <- .initializeEmptySppList(adduct, spwin)


                    #// de-protonated ----
                    if("M-H" %in% adduct) {

                          lipTmp    <- sldb$`Exact m/z of [M-H]-`[i]

                          if(!is.na(lipTmp)) { # for example there is no protonated version of this lipid species

                                # metaData list
                                mt  <- list(mode = "negative",
                                            adduct = "M-H",
                                            lipidID = sldb$`Lipid ID`[i],
                                            sumformula = sldb$`Formula (pH7.3)`[i],
                                            abbrev = sldb$`Abbreviation*`[i],
                                            numDoubleBonds = sldb$numDoubleBond[i],
                                            lipidClass = sldb$lipidGroup[i],
                                            chainLength = sldb$chainLength[i])

                                sppCotainer[["M-H"]] <- searchAnalyte(m = lipTmp,
                                                          fwhm = getFwhm(fwhmObj, lipTmp),
                                                          spData = spData,
                                                          wMethod = wMethod,
                                                          verifiedMasses = verifiedMasses,
                                                          confirmedOnly = confirmedOnly,
                                                          metaData = mt)

                                #sppCotainer <- c(sppCotainer, list(hitsDeprot))

                                rm(mt, lipTmp)


                          }

                    }


                    #// protonated ----
                    if("M+H" %in% adduct){
                          lipTmp        = sldb$`Exact m/z of [M+H]+`[i]

                          if(!is.na(lipTmp)) { # for example there is no protonated version of this lipid species


                                # metaData list
                                mt  <- list(mode = "positive",
                                            adduct = "M+H",
                                            lipidID = sldb$`Lipid ID`[i],
                                            sumformula = sldb$`Formula (pH7.3)`[i],
                                            abbrev = sldb$`Abbreviation*`[i],
                                            numDoubleBonds = sldb$numDoubleBond[i],
                                            lipidClass = sldb$lipidGroup[i],
                                            chainLength = sldb$chainLength[i])

                                sppCotainer[["M+H"]] <- searchAnalyte(m = lipTmp,
                                                        fwhm = getFwhm(fwhmObj, lipTmp),
                                                        spData = spData,
                                                        wMethod = wMethod,
                                                        verifiedMasses = verifiedMasses,
                                                        confirmedOnly = confirmedOnly,
                                                        metaData = mt)

                                #sppCotainer <- c(sppCotainer, list(hitsProt))


                                rm(mt, lipTmp)



                          }
                    }



                    #// Na+ adduct ----
                    if("M+Na" %in% adduct){
                          lipTmp               = sldb$`Exact m/z of [M+Na]+`[i]

                          if(!is.na(lipTmp)) { # for example there is no Na-adduct version



                                # metaData list
                                mt  <- list(mode = "positive",
                                            adduct = "M+Na",
                                            lipidID = sldb$`Lipid ID`[i],
                                            sumformula = sldb$`Formula (pH7.3)`[i],
                                            abbrev = sldb$`Abbreviation*`[i],
                                            numDoubleBonds = sldb$numDoubleBond[i],
                                            lipidClass = sldb$lipidGroup[i],
                                            chainLength = sldb$chainLength[i])

                                sppCotainer[["M+Na"]] <- searchAnalyte(m = lipTmp,
                                                       fwhm = getFwhm(fwhmObj, lipTmp),
                                                       spData = spData,
                                                       wMethod = wMethod,
                                                       verifiedMasses = verifiedMasses,
                                                       confirmedOnly = confirmedOnly,
                                                       metaData = mt)

                                #sppCotainer <- c(sppCotainer, list(hitsSod))


                                rm(mt, lipTmp)



                          }

                    }


                    #// K+ adduct ----
                    if("M+K" %in% adduct | "M+k" %in% adduct){

                          lipTmp        = sldb$`Exact m/z of [M+K]+`[i]

                          if(!is.na(lipTmp)) { # for example there is no Na-adduct version

                                # metaData list
                                mt  <- list(mode = "positive",
                                            adduct = "M+K",
                                            lipidID = sldb$`Lipid ID`[i],
                                            sumformula = sldb$`Formula (pH7.3)`[i],
                                            abbrev = sldb$`Abbreviation*`[i],
                                            numDoubleBonds = sldb$numDoubleBond[i],
                                            lipidClass = sldb$lipidGroup[i],
                                            chainLength = sldb$chainLength[i])

                                sppCotainer[["M+K"]] <- searchAnalyte(m = lipTmp,
                                                         fwhm = getFwhm(fwhmObj, lipTmp),
                                                         spData = spData,
                                                         wMethod = wMethod,
                                                         verifiedMasses = verifiedMasses,
                                                         confirmedOnly = confirmedOnly,
                                                         metaData = mt)

                               #sppCotainer <- c(sppCotainer, list(hitsPotas))


                                rm(mt, lipTmp)



                          }
                    }

                    # noDetections <- sapply(sppCotainer, is.null)
                    #
                    # if(all(!noDetections)){ # if there are no detections return empty ppp object
                    #       return(spatstat.geom::ppp(x = integer(0), y = integer(0), window = spwin))
                    # } else{
                    #       return(superimposeAnalytes(sppCotainer, spWin = spwin))
                    # }

                    return(superimposeAnalytes(sppCotainer, spWin = spwin))

                  })




            #// merge
            hitsList <- superimposeAnalytes(hitsList, spWin = spwin)

            # if(confirmedOnly){ # manually set this flag
            #       hitsList$metaData$mzConfirmed <- TRUE
            # }


      # if(verbose)
      #       close(pb)


      return(hitsList)


}


# .reduceSearchList <- function(searchList){ # remove classes which did not generate hits
#
#       tokeep <- !(sapply(searchList$hitsList, function(x) spatstat.geom::is.empty.ppp(x)))
#
#       searchList <- lipidSearchList(lipidList = searchList$lipidList,
#                                     hitsList = searchList$hitsList[tokeep],
#                                     allClasses = names(searchList$hitsList)[tokeep])
#
#       return(searchList)
#
# }

.initializeEmptySppList <- function(vec, spwin){

      # this is used to initialize an empty list of spp objects
      # vec is a character vector holding the query adducts
      # spwin is the spatial window of type spatstat.geom::owin

      emptyspp <- spatstat.geom::ppp(x = integer(0), y = integer(0),
                                     marks = data.frame(idx = integer(0), intensity = numeric(0)),
                                     window = spwin, checkdup = FALSE, drop = FALSE)

      cont <- replicate(length(vec), emptyspp, simplify = FALSE)

      names(cont) <- vec

      return(cont)

}

.filtersldb <- function(sldb, verifiedMasses, adduct, fwhmObj){

      # filter sldb to include only verified masses - speed up computations

      # to fix: consider adding more stringent filtration on verifiedMasses to also
      # includ the adduct from which these masses were calculated.

      idxToKeep <- integer(0)

      if("M-H" %in% adduct){

            tol <- getFwhm(fwhmObj, sldb$`Exact m/z of [M-H]-`)

            if(any(is.na(tol))){
                  tol[is.na(tol)] <- 0
            }


            tmp <- MALDIquant::match.closest(x = sldb$`Exact m/z of [M-H]-`,
                                              table = sort(verifiedMasses),
                                              tolerance = tol)
            idxToKeep <- c(idxToKeep, which(!is.na(tmp)))
      }

      if("M+H" %in% adduct){

            tol <- getFwhm(fwhmObj, sldb$`Exact m/z of [M+H]+`)

            if(any(is.na(tol))){
                  tol[is.na(tol)] <- 0
            }

            tmp <- MALDIquant::match.closest(x = sldb$`Exact m/z of [M+H]+`,
                                              table = sort(verifiedMasses),
                                              tolerance = tol)
            idxToKeep <- c(idxToKeep, which(!is.na(tmp)))

      }

      if("M+Na" %in% adduct){

            tol <- getFwhm(fwhmObj, sldb$`Exact m/z of [M+Na]+`)

            if(any(is.na(tol))){
                  tol[is.na(tol)] <- 0
            }

            tmp <- MALDIquant::match.closest(x = sldb$`Exact m/z of [M+Na]+`,
                                             table = sort(verifiedMasses),
                                             tolerance = tol)
            idxToKeep <- c(idxToKeep, which(!is.na(tmp)))
      }

      if("M+K" %in% adduct){

            tol <- getFwhm(fwhmObj, sldb$`Exact m/z of [M+K]+`)

            if(any(is.na(tol))){
                  tol[is.na(tol)] <- 0
            }

            tmp <- MALDIquant::match.closest(x = sldb$`Exact m/z of [M+K]+`,
                                             table = sort(verifiedMasses),
                                             tolerance = tol)
            idxToKeep <- c(idxToKeep, which(!is.na(tmp)))
      }

      if(length(idxToKeep) == 0) {
            stop("No matches of 'verifiedMasses' in the SwissLipid DB. \n" )
      }

      sldb <- sldb[idxToKeep, ]


}

