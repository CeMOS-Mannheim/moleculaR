#' Load The Swisslipids database
#'
#' This function reads the swisslipids database previously douwnloaded as a file. It keeps the entries at "Species" level and
#' adds a columns for the number of double bonds. 
#'
#' @param path: 	path to the swiss database file.
#' @return
#' a filtered data frame holding the entries of the swisslipids database at the "species" level. 
#'
#' @export
#'
loadSwissDB    = function(path) {


		##// load swisslipids database ----
		swiss                = read.csv(file = path, 
		                                   sep = '\t', header = TRUE, stringsAsFactors = F,
		                                 encoding = "UTF-8", check.names = F)
		 
		 
		# focus on species because we will not be able to distinguish between subspecies anyway 
		swiss                = swiss[swiss$Level == "Species", ]
 
 
		# remove the entries with no abbreviation
		swiss                = swiss[swiss$`Abbreviation*` != "", ]
 
 
		# create a column for the number of double bonds (C=C)
		swiss$numDoubleBond  = vector("integer", nrow(swiss))
 
		tmp                  = sub(".*:", "", swiss$`Abbreviation*`)
		swiss$numDoubleBond  = as.integer(sub(")", "", tmp))
 
		rm(tmp)


       

       return(swiss)

  }

