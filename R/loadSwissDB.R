#' Load The Swisslipids database
#'
#' This function reads the swisslipids database previously douwnloaded as a file. It keeps the entries at "Species" level and
#' adds a columns for the number of double bonds.
#'
#' @param path: 	path to the SwissLipids database file.
#' @return
#' a filtered data frame holding the entries of the swisslipids database at the "species" level. Entries with no
#' `Abbreviation*` record are removed. Two additional columns are added `numDoubleBond` and `lipidClass` to hold
#' the calculated number of double-bonds and the lipid class, respectively. This is done to make easier to
#' filter for these attriutes.
#'
#' @export
#'
loadSwissDB    <- function(path) {


		##// load swisslipids database ----
		sldb                <- read.csv(file = path,
		                                   sep = '\t', header = TRUE, stringsAsFactors = F,
		                                 encoding = "UTF-8", check.names = F)


		# focus on species because we will not be able to distinguish between subspecies anyway
		sldb                <- sldb[sldb$Level == "Species", ]


		# remove the entries with no abbreviation
		sldb                <- sldb[sldb$`Abbreviation*` != "", ]


		# add number of double bonds info as an additional column to sldb
		sldb                <- .appendNumDoubleBond(sldb)


		# add lipid class info as an additional column to sldb
		sldb                <- .appendLipidClass(sldb)

		# add fatty acid chain length as an additional column to sldb
		sldb                <- .appendChainLength(sldb)



       return(sldb)

}

.appendNumDoubleBond <- function(sldb){

      # This appends the number of double bonds info as an additional column
      # to the SwissLipids database `sldb`

      numDoubleBond       <- vector("integer", nrow(sldb))

      tmp                 <- sub(".*:", "", sldb$`Abbreviation*`)

      numDoubleBond       <- as.integer(sub(")", "", tmp))

      return(cbind(sldb, numDoubleBond = numDoubleBond))

}

.appendLipidClass <- function(sldb){

      # This appends the lipid class info in the form of PC(x:x) as an additional column
      # to the SwissLipids database `sldb`

      allAbbrevs           <- strsplit(x = sldb$`Abbreviation*`, split = "\\(")
      allAbbrevs           <- sapply(allAbbrevs, FUN = function(i) paste0(c(i[1], "(", gsub("[[:digit:]]", "", i[[2]])), collapse = ""))
      allAbbrevs           <- strsplit(allAbbrevs, ":")
      allAbbrevs           <- sapply(allAbbrevs, FUN = function(i) paste0(i[1], "x:x" , i[2], collapse = ""))

      return(cbind(sldb, lipidGroup = allAbbrevs))
}

.appendChainLength <- function(sldb) {

      # This appends the length of the fatty acid chain as an additional column
      # to the SwissLipids database `sldb`

      allAbbrevs        <- strsplit(x = sldb$`Abbreviation*`, split = "\\(")
      allAbbrevs        <- sapply(allAbbrevs, FUN = function(i) strsplit(x=i[2], ":"))
      allAbbrevs        <- as.integer(sapply(allAbbrevs, FUN = function(i) gsub(pattern = "[^0-9]", replacement = "", x = i[1])))

      return(cbind(sldb, chainLength = allAbbrevs))
}

