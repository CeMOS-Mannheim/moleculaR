#' Fix metaspace entries
#'
#' This function fixes metaspace entires according to envipat package. 
#'
#' @param df: 	the dataframe holding the metaspace identification data
#' @return
#' a fixed dataframe.  
#'
#' @export
#'
fixMetaSpace    = function(df) {
       
       # fix sumforulas and adducts representation of metaspace - make easier to compare with hits
       data("isotopes", package = "enviPat")
       data("adducts", package = "enviPat")

       df$formula = enviPat::check_chemform(isotopes = isotopes,
                                            chemforms = df$formula,
                                            get_sorted = TRUE)$new_formula


       df$adduct = gsub(x = df$adduct,
                        pattern = "M-",
                        replacement = "",
                        ignore.case = FALSE)
       
       df$adduct = enviPat::check_chemform(isotopes = isotopes,
                                           chemforms = df$adduct)$new_formula

       
       
       
       
       
       return(df)
       
}



