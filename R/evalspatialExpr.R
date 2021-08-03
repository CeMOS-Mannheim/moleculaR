#' Evaluates an expression of SPP objects
#'
#' This evaluates an arbitrary expression involving spatial point pattern objects 'spatstat::ppp' which encode
#' analyte distributions of MSI data.
#'
#' @param expr: 	   An unquoted (unwrapped) expression of 'ppp' objects listed in 'dataList'.
#' @param dataList:  A named list containing the 'ppp' objects involved in 'expr'.
#' @param ppwin:     The window 'spatstat::owin' object which all elements of 'dataList' share.
#'
#' @return
#' A list containing 'spatstat::im' and the corresponding 'spatstat::ppp' objects resulting of
#' the applied expression 'expr'.
#'
#' @export
#'

evalSpatialExpr <- function(exprn, dataList, ppwin, bw, sqrtTransform = FALSE){

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
            dataList[[ivar]] <- spatstat::pixellate(dataList[[ivar]],
                                      weights = dataList[[ivar]]$marks$intensity,
                                      W = spatstat::as.mask(ppwin,
                                                            dimyx=c(diff(ppwin$yrange),
                                                                    diff(ppwin$xrange))),
                                      padzero = FALSE, savemap = FALSE)


      }

      # the result of the expression as an 'im' object
      res <- eval(expr = e, envir = dataList)
      #res <- spatstat::eval.im(exprn, envir = dataList)

      # convert back to ppp
      imdf <- spatstat::as.data.frame.im(x = res)

      if(any(is.infinite(imdf$value))){
            imdf <- imdf[which(!is.infinite(imdf$value)), ]
      }

      exprppp <- spatstat::ppp(x = imdf$x, y = imdf$y,
                               marks = data.frame(intensity = imdf$value),
                               window = ppwin, checkdup = FALSE, drop = FALSE)

      exprim <- spatstat::pixellate(exprppp,
                                    weights = exprppp$marks$intensity,
                                    W = spatstat::as.mask(ppwin,
                                                          dimyx=c(diff(ppwin$yrange) + 0,
                                                                  diff(ppwin$xrange) + 0)),
                                    padzero = FALSE, savemap = FALSE)

      return(list(exprim = exprim, exprppp = exprppp))




}
