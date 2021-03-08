#' Create owin window lists
#'
#' This function creates `spatstat::owin` lists of polygonal information stored in the mis file (regions). 
#' This function only focuses on one tissue region and the corresponding annotation regions.
#'
#' @param regions: 	The regions object as created by `msiTools::GetRegions`.
#' @param regName: The indices pointing to the tissue regions. 
#' @return
#' A list of i) tissue window objects and ii) annotations window objects.   
#'
#' @export
#'
createWinList2        = function(regions, regName, plotIt = FALSE ) {
       
       
       #// checks
       atiss         = spatstat.utils::Area.xypolygon(regions$polyPoints[[regName]])
       #a             = sapply(regions$polyPoints[[regName]], spatstat.utils::Area.xypolygon)
       
       if(atiss < 0)
              stop("The given tissue polygon is traversed in the wrong direction - according to spatstat definition.\n")
              
       
       #// fix the regions representation into x,y and area - done to make idintifying holes (indicated by -ve area) easier 
       fixedRegions  = lapply(regions$polyPoints, FUN = function(i) {
              
              if(nrow(i) < 3) {
                     
                     return(list(x = i$x, y = i$y, area = abs(diff(i$x) * diff(i$y))))
              }
              
              return(list(x = i$x, y = i$y, area = spatstat.utils::Area.xypolygon(i)))
       })
       
              
       #// create win objects from regions of interest
       winTissue      =  spatstat::owin(poly = fixedRegions[[regName]])
       
       
       #// creat owin objects from the annotations in the negative ion mode mis file
       #   note that the annotations are stored in the negative mis file i.e. regionsNeg
       
       winAnnot      = setNames(vector("list", 5), c("MVP", "LE", "nec", "CT", "IT"))
       
       for(ilist in names(winAnnot)) {
              
              idx    = grep(pattern = ilist, 
                            x = regions$minMaxCoords$ROI_name, ignore.case = FALSE)
              
              if(length(idx) == 0) {
                     winAnnot[[ilist]] = NULL
                     next
              }
              
              # add a for loop to loop through all annotation regions
              
              #// put limits on the x and y directions
              # create a rectangular owin object and test against it
              rectWin       = spatstat::owin(xrange=winTissue$xrange + c(-20,20), 
                                             yrange=winTissue$yrange + c(-20,20))
              
              inWin         = sapply(idx, FUN = function(i) {
                     
                     any(spatstat::inside.owin(x = fixedRegions[[i]]$x, 
                                           y = fixedRegions[[i]]$y, 
                                           w = rectWin))
              })
              

              if(any(inWin)) {
                    
                     
                     winAnnot[[ilist]]     = spatstat::owin(poly = fixedRegions[idx[inWin]]) 
                     
              }
              
              
              
              
       }
       
       if(plotIt) {
              
              cols          = c(MVP = "red", LE = "yellow", nec = "black", CT = "blue", IT = "green" )
              
              spatstat::plot.owin(winTissue, ylim = rev(range(winTissue$y)))
              
              foo           = sapply(seq(1,length(winAnnot)), FUN = function(i) {
                     
                     if(!is.null(winAnnot[[i]])) {
                            
                            spatstat::plot.owin(winAnnot[[i]], 
                                                add = TRUE, col = cols[i])
                            
                     }
                            
              })
              
              
       }
       
       
       
       return(list(winTissue = winTissue, 
                   winAnnot = winAnnot))
       
       
       
       
}


# spPolyList           = lapply(annotateIdx, FUN = function(x) { #import via sp package
#        
#        sp::Polygon(coords = e$fixedRegions[[x]])
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


