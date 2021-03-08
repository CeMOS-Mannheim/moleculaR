#' Calculates isotop patterns for a given sumformula
#'
#' 
#'
#' @param sumFormula:       not much to say here. 
#' @param charge:           \code{z} in \code{m/z}.
#' @param adduct:           the expected adduct formation.  
#' @param subtracted:       if the adduct is added (ex. protonated) or subtracted (ex. substracted).
#' @param plot:             whether to plot the resulting isotop patterns. : 
#' @return
#' A dataframe containing the calculated isotop patterns. 
#'
#' @export
#'
sumFormToMass        = function(sumFormula, charge = 1L, adduct = "H1", subtracted = FALSE, plot = FALSE) {
       
       
       
       data("isotopes", package = "enviPat")
       data("adducts", package = "enviPat")
       
       sumFormula    = enviPat::check_chemform(isotopes = isotopes, 
                                               chemforms = sumFormula,
                                               get_sorted = FALSE)$new_formula
       
       
       
       
       #// a function to compute the m/z for the sumformula for each analyte using EnviPat. 
       
       if(subtracted) 
       {
               txt    = paste0(sumFormula, " - ", adduct)
              sumFormula = enviPat::subform(sumFormula, adduct)
              
              
       } else
       {
               txt    = paste0(sumFormula, " + ", adduct)
              sumFormula = enviPat::mergeform(sumFormula, adduct)
              
              
       }
       
       patterns      = enviPat::isopattern(isotopes = isotopes, chemforms = sumFormula, threshold = 1,
                                           charge = charge, plotit = FALSE, verbose = FALSE)[[1]]
       
       if(plot) {
              
              plot(x = patterns[ , 1], y = patterns[ , 2], type = "h", bty = "n", xaxt = "n", 
                   main = txt, xlab = "m/z", ylab = "relative abundance")
              axis(side = 1, at = patterns[,1])
              
       }
       
              
              
       return(patterns)
       
       
       
       
}




# #// a function to compute the m/z for the sumformula for each analyte using EnviPat. 
# isopat        = function(form, charge, adduct, subtracted = FALSE, plot = FALSE)
# {
#        
#        #data("isotopes", package = "enviPat")
#        #data("adducts", package = "enviPat")
#        
#        #form          = enviPat::check_chemform(isotopes, form)
#        #form          = enviPat::multiform(form[ , 2],1)
#        
#        if(subtracted) 
#        {
#               form   = enviPat::subform(form, adduct)
#               
#        } else
#        {
#               form   = enviPat::mergeform(form, adduct)
#               
#        }
#        
#        patterns      = enviPat::isopattern(isotopes, form, charge = charge, plotit = plot, verbose = FALSE)
#        
#        return(patterns)
#        
# }

