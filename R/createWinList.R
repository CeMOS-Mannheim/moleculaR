#' Create owin window lists
#'
#' This function creates `spatstat::owin` lists of polygonal information stored in the mis file (regions). 
#'
#' @param regions: 	The regions object as created by `msiTools::GetRegions`.
#' @param regIdx: The indices pointing to the tissue regions. 
#' @return
#' A list of i) tissue window objects and ii) annotations window objects.   
#'
#' @export
#'
createWinList        = function(regions, regIdx ) {
       
       #// create win objects from regions of interest
       winList       = setNames(object = lapply(regIdx, FUN = function(x) spatstat::owin(poly = regions$polyPoints[[x]])), 
                                   nm = names(regions$polyPoints)[regIdx])
       
       
       #// creat owin objects from the annotations in the negative ion mode mis file
       #   note that the annotations are stored in the negative mis file i.e. regionsNeg
       annotationList= setNames(vector("list", 5), c("MVP", "LE", "nec", "CT", "IT"))
       
       for(ilist in names(annotationList)) {
              
              idx    = grep(pattern = ilist, 
                            x = regions$minMaxCoords$ROI_name, ignore.case = FALSE)
              
              if(length(idx) == 0) {
                     annotationList[[ilist]] = NULL
                     next
              }
              
              annotationList[[ilist]]     = spatstat::owin(poly = regions$polyPoints[idx])
              
              
       }
       
       
       
       return(list(tissueWin = winList, 
                   annotationWin = annotationList))
       
       
       

}


# spPolyList           = lapply(annotateIdx, FUN = function(x) { #import via sp package
#        
#        sp::Polygon(coords = e$regions$polyPoints[[x]])
# })
# annotationsWin       = spatstat::owin(poly = e$regionsNeg$polyPoints[annotateIdx])
# 
# annotationsList      =  lapply(annotateIdx, FUN = function(i) { 
#        
#        print(i)
#        x             = e$regionsNeg$polyPoints[[i]]$x
#        y             = e$regionsNeg$polyPoints[[i]]$y
#        isCCW         = ifelse(sum(diff(x) * (head(y, -1) + tail(y, -1))) < 0, TRUE, FALSE)
#        
#        if(isCCW) {return(NULL)}
#         
#        spatstat::owin(p


