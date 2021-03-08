bw.scott <- function(X, isotropic=FALSE, d=NULL) {
       #stopifnot(is.ppp(X) || is.lpp(X) || is.pp3(X) || is.ppx(X))
       if(is.null(d)) { d <- spatstat::spatdim(X) } else check.1.integer(d)
       nX <- spatstat::npoints(X)
       cX <- spatstat::coords(X, spatial=TRUE, temporal=FALSE, local=FALSE)
       sdX <- apply(cX, 2, sd)
       if(isotropic) {
              #' geometric mean
              sdX <- exp(mean(log(pmax(sdX, .Machine$double.eps))))
       }
       b <- sdX * nX^(-1/(d+4))
       names(b) <- if(isotropic) "sigma" else paste0("sigma.", colnames(cX))
       return(b)
}

bw.scott.iso <- function(X) { bw.scott(X, isotropic=TRUE) }