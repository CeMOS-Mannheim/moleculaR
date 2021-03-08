#' Add annotation contours to a plot
#'
#' This function creates an analyte probability map from a given spatial point pattern. 
#'
#' @param winList:      A list of owin objects containing the annotation contours. 
#' @param .alpha:    Alpha as a percentage
#' @param lwd:          line width
#' @param lty:          line type
#' @return
#' Nothing. plots lines.  
#'
#' @export
#'
addAnnotations          = function(winList, .alpha = 0.5, lwd = 1, lty = "solid") {
      
      
   
   
      # col               = as.vector(col2rgb("blue"))
      # blue              = rgb(col[1], col[2], col[3], alpha = alpha, max = 255)
      # 
      # col               = as.vector(col2rgb("green"))
      # green             = rgb(col[1], col[2], col[3], alpha = alpha, max = 255)
      # 
      # col               = as.vector(col2rgb("red"))
      # red               = rgb(col[1], col[2], col[3], alpha = alpha, max = 255)
      # 
      # col               = as.vector(col2rgb("yellow"))
      # yellow            = rgb(col[1], col[2], col[3], alpha = alpha, max = 255)
      # 
      # col               = as.vector(col2rgb("black"))
      # black             = rgb(col[1], col[2], col[3], alpha = alpha, max = 255)
      # 
      
      for(ix in names(winList)) {
            
            if(is.null(winList[[ix]])) {next}
            
            xcol    = switch(ix, 
                             CT = t_col("blue", .alpha = .alpha),
                             IT = t_col("green", .alpha = .alpha), 
                             MVP = t_col("red", .alpha = .alpha), 
                             LE = t_col("yellow", .alpha = .alpha), 
                             nec = t_col("black", .alpha = .alpha))
            
            spatstat::plot.owin(winList[[ix]], 
                                border = xcol, 
                                lty = lty, 
                                lwd = lwd, 
                                add = TRUE) 
            
      }
      
     
}



## Transparent colors

t_col <- function(color, .alpha = 0.5) {

   
   ## Get RGB values for named color
   rgb.val <- col2rgb(color)
   
   ## Make new color using input color as base and alpha set by transparency
   t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
                max = 255,
                alpha = .alpha * 255)
   
   ## Save the color
   invisible(t.col)
}
## END