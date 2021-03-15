#' Initialize lipids search
#'
#' This function initializes some lists for lipid search.
#'
#' @param swissdb: 	The swiss lipids db as dataframe.
#' @return
#' a list of two entreis; a list containing segmented lipid
#' species into individual entries and an additional empty list
#' with the same structure and naming to store the search results.
#'
#' @export
#'
initLipidSearch     = function(swissdb) {



       allAbbrevs           = strsplit(x = swissdb$`Abbreviation*`, split = "\\(")
       allAbbrevs           = sapply(allAbbrevs, FUN = function(i) paste0(c(i[1], "(", gsub("[[:digit:]]", "", i[[2]])), collapse = ""))
       allAbbrevs           = strsplit(allAbbrevs, ":")
       allAbbrevs           = sapply(allAbbrevs, FUN = function(i) paste0(i[1], "x:x" , i[2], collapse = ""))


       allGroups            = unique(allAbbrevs)
       swissList            = setNames(object = vector("list", length = length(allGroups)), nm = allGroups)
       lipidHits            = setNames(object = vector("list", length = length(allGroups)), nm = allGroups)

       # #// focus on the alphanumeric part of the abbreviation
       # allAbbrevs           = gsub(pattern = "[^[:alpha:]]", replacement = "", x = swissdb$`Abbreviation*`)
       #


       for(i in allGroups)
       {

              #idx           = grep(pattern = paste0("^", i, "$"), x = allAbbrevs, ignore.case = FALSE)
              idx           = grep(pattern = i, x = allAbbrevs, fixed = TRUE)

              swissList[[i]]= swissdb[idx, ]

       }

       # O- -> ethers || nothing of a- -> acyl groups or esters



       return(list(swissList = swissList,
                   lipidHits = lipidHits,
                   allClasses = allGroups))

}

