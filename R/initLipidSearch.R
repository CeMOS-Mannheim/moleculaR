#' Initialize lipids search
#'
#' This function initializes some lists for lipid search.
#'
#' @param swissdb: 	The swiss lipids db as dataframe.
#' @return
#' an S3 object with three entreis; a list containing lipid entries
#' organized into lipid classes, an additional empty 'hitslist'
#' with the same structure and naming to store the search results and
#' a character verctor summarizing all lipid classes in the first entry.
#'
#' @export
#'
initLipidSearch     <- function(swissdb) {



       allAbbrevs           <- strsplit(x = swissdb$`Abbreviation*`, split = "\\(")
       allAbbrevs           <- sapply(allAbbrevs, FUN = function(i) paste0(c(i[1], "(", gsub("[[:digit:]]", "", i[[2]])), collapse = ""))
       allAbbrevs           <- strsplit(allAbbrevs, ":")
       allAbbrevs           <- sapply(allAbbrevs, FUN = function(i) paste0(i[1], "x:x" , i[2], collapse = ""))


       allGroups            <- unique(allAbbrevs)
       swissList            <- setNames(object = vector("list", length = length(allGroups)), nm = allGroups)
       lipidHits            <- setNames(object = vector("list", length = length(allGroups)), nm = allGroups)



       allAbbrevsUnpunct   <- gsub("\\(|\\)", " ", allAbbrevs)

       for(i in allGroups)
       {

              ix            <- gsub("\\(|\\)", " ", i)
              idx           <- grep(pattern = paste0("^", ix, "$"), x = allAbbrevsUnpunct)

              swissList[[i]]<- swissdb[idx, ]

       }



       return(lipidSearchList(lipidList = swissList,
                              hitsList = lipidHits,
                              allClasses = allGroups))
}

