#' Evaluates an expression of SPP objects
#'
#' This evaluates an arbitrary expression involving spatial point pattern objects 'ppp' which encode
#' analyte distributions of MSI data.
#'
#' @param expr: 	   A character vector representing tha arithmetic expression of 'ppp' objects listed in 'dataList'.
#' @param dataList:  A named list containing the 'ppp' objects involved in 'expr'.
#' @param ppwin:     The window 'owin' object which all elements of 'dataList' share.
#' @param sqrtTransform: Optional square root transofmation of the 'ppp' objects in 'dataList'.
#'
#' @return
#' An object of type 'ppp' and 'moleculaR::analytePointPattern' holding the result of
#' the applied expression 'expr'.
#'
#' @export
#' @include manualSpatstatImport.R

evalSpatialExpr <- function(exprn, dataList, ppwin, sqrtTransform = FALSE){

      e <- parse(text = exprn)

      vars <- all.vars(e)

      if(sqrtTransform){
            dataList <- lapply(dataList, function(i) {
                  i$marks$intensity <- sqrt(i$marks$intensity)
                  return(i)
            })
      }


      for(ivar in vars){
            if(!(ivar %in% names(dataList)))
                  stop(paste0(ivar, " is not found in 'dataList'. \n"))


            # create a pixellated image and assign back
            dataList[[ivar]] <-  spp2im(dataList[[ivar]])


      }

      # the result of the expression as an 'im' object
      res <- eval(expr = e, envir = dataList)
      #res <- eval.im(exprn, envir = dataList)

      # convert back to ppp
      imdf <- as.data.frame.im(x = res)

      if(any(is.infinite(imdf$value))){
            imdf <- imdf[which(is.finite(imdf$value)), ]
      }


      exprppp <- analytePointPattern(x = imdf$x, y = imdf$y,
                                     intensity = imdf$value,
                                     win = ppwin, mzVals = "expr")



      return(exprppp)




}
