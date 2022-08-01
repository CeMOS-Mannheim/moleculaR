#' Evaluates an expression of SPP objects
#'
#' This evaluates an arbitrary expression involving spatial point pattern objects 'ppp' which encode
#' analyte distributions of MSI data.
#'
#' @param expr: 	   A character vector representing the arithmetic expression of 'ppp' objects listed in 'dataList'.
#' @param dataList:  A named list containing the 'ppp' objects involved in 'expr'.
#' @param win:     The window 'owin' object which all elements of 'dataList' share.
#' @param normMethod: a character, the method used for normalizing intensities of `ppp` objects
#' in `dataList`. This is particularly important when `CPPMs` are considered i.e. `ppp` objects
#' carrying more than one analyte. One of `c("imageMean", "pixelMean")`. Defaults to `imageMean`.
#'
#' @details
#' For `normMethod='imageMean'` if the terms of `expr` are CPPMs (i.e. holding
#' several analytes) each pixel intensity within CPPM image is divided by the
#' number of analytes contributing to the CPPM. For `normMethod='pixelWise'`,
#' each pixel is down-weighted by the absolute number of analytes present in it.
#' Both normalization types are applied before applying the arithmetic expression
#' and do not have any effect of `spp`s that hold a single analyte.
#'
#' @return
#' An object of type 'ppp' and 'moleculaR::analytePointPattern' holding the result of
#' the applied expression 'expr'.
#'
#' @export
#' @include manualSpatstatImport.R

evalSpatialExpr <- function(exprn, dataList, win, normMethod = "imageMean"){

      e <- parse(text = exprn)

      vars <- all.vars(e)


      for(ivar in vars){
            if(!(ivar %in% names(dataList)))
                  stop(paste0(ivar, " is not found in 'dataList'. \n"))


            dataList[[ivar]] <- switch(normMethod,
                                       "imageMean" = {
                                             # Normalize intensities of each image by the number
                                             # of mz values (analytes) that compose that image.
                                             spp2im(dataList[[ivar]]) / length(dataList[[ivar]]$metaData$mzVals)
                                       },
                                       "pixelMean" = {
                                             # Normalize intensities of each pixel in an image by the number
                                             # of mz values (analytes) that contributes to that pixel.
                                             px <- spp2im(dataList[[ivar]], weighted = FALSE) # stores pixel counts
                                             spp2im(dataList[[ivar]]) / px
                                       },
                                       stop("incorrect input for 'normMethod'.\n"))



      }

      # the result of the expression as an 'im' object
      res <- eval(expr = e, envir = dataList)
      #res <- eval.im(exprn, envir = dataList)

      # sanity checks
      imdf <- as.data.frame.im(x = res)

      if(any(is.infinite(imdf$value))){
            imdf <- imdf[which(is.finite(imdf$value)), ]
      }

      if(any(imdf$value == 0)){
            imdf <- imdf[which(imdf$value != 0), ]
      }

      # convert back to spp
      exprppp <- analytePointPattern(x = imdf$x, y = imdf$y,
                                     intensity = imdf$value,
                                     win = win, mzVals = "expr")



      return(exprppp)




}
