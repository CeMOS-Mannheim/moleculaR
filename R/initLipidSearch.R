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
       
       
       
       swissList            = list(PC = NULL, 
                                   LPC = NULL, 
                                   PE = NULL, 
                                   LPE = NULL, 
                                   PI = NULL, 
                                   LPI = NULL,
                                   PS = NULL, 
                                   LPS = NULL,
                                   TG = NULL, 
                                   SM = NULL,
                                   PG = NULL, 
                                   PA = NULL, 
                                   DG = NULL, 
                                   PIP = NULL, 
                                   PIP2 = NULL, 
                                   PIP3 = NULL ,
                                   FA = NULL)
       
       lipidHits             = list(PC = NULL, 
                                    LPC = NULL, 
                                    PE = NULL, 
                                    LPE = NULL, 
                                    PI = NULL, 
                                    LPI = NULL,
                                    PS = NULL, 
                                    LPS = NULL,
                                    TG = NULL, 
                                    SM = NULL,
                                    PG = NULL, 
                                    PA = NULL, 
                                    DG = NULL, 
                                    PIP = NULL, 
                                    PIP2 = NULL, 
                                    PIP3 = NULL,
                                    FA = NULL)
       
       
       
       # O- -> ethers || nothing of a- -> acyl groups or esters
       
       #// focus on the alphanumeric part of the abbreviation
       #allAbbrevs           = gsub(pattern = "[^[:alpha:]]", replacement = "", x = swissdb$`Abbreviation*`)
       allAbbrevs           = strsplit(x = swissdb$`Abbreviation*`, split = "\\(")
       allAbbrevs           = sapply(allAbbrevs, FUN = function(i) paste0(c(i[1], gsub("[^[:alpha:]]", "", i[[2]])), collapse = "-"))
       
       # filter for PC
       idx                  = grep(pattern = "^PC-$", x = allAbbrevs, ignore.case = FALSE)
       swissList$PC         = swissdb[idx, ]
       
       # filter for LPC
       idx                  = grep(pattern = "^LPC-$", x = allAbbrevs, ignore.case = FALSE)
       swissList$LPC        = swissdb[idx, ]
       
       
       # filter for PE
       idx                  = grep(pattern = "^PE-$", x = allAbbrevs, ignore.case = FALSE)
       swissList$PE        = swissdb[idx, ]
       
       # filter for LPE
       idx                  = grep(pattern = "^LPE-$", x = allAbbrevs, ignore.case = FALSE)
       swissList$LPE        = swissdb[idx, ]
       
       
       # filter for PI
       idx                  = grep(pattern = "^PI-$", x = allAbbrevs, ignore.case = FALSE)
       swissList$PI         = swissdb[idx, ]
       
       # filter for LPI
       idx                  = grep(pattern = "^LPI-$", x = allAbbrevs, ignore.case = FALSE)
       swissList$LPI        = swissdb[idx, ]
       
       # filter for PS
       idx                  = grep(pattern = "^PS-$", x = allAbbrevs, ignore.case = FALSE)
       swissList$PS         = swissdb[idx, ]
       
       # filter for LPS
       idx                  = grep(pattern = "^LPS-$", x = allAbbrevs, ignore.case = FALSE)
       swissList$LPS        = swissdb[idx, ]
       
       
       # filter for TG
       idx                  = grep(pattern = "^TG-$", x = allAbbrevs, ignore.case = FALSE)
       swissList$TG         = swissdb[idx, ]
       
       # filter for SM
       idx                  = grep(pattern = "^SM-", x = allAbbrevs, ignore.case = FALSE)
       swissList$SM         = swissdb[idx, ]
       
       # filter for PG
       idx                  = grep(pattern = "^PG-$", x = allAbbrevs, ignore.case = FALSE)
       swissList$PG         = swissdb[idx, ]
       
       # filter for PA
       idx                  = grep(pattern = "^PA-$", x = allAbbrevs, ignore.case = FALSE)
       swissList$PA         = swissdb[idx, ]
       
       # filter for DG
       idx                  = grep(pattern = "^DG-$", x = allAbbrevs, ignore.case = FALSE)
       swissList$DG         = swissdb[idx, ]
       
       # filter for PIP
       # idx                  = grep(pattern = "^(PIP)", x = swissdb$`Abbreviation*`, ignore.case = FALSE)
       # swissList$PIP        = swissdb[idx, ]
       # idx                  = grep(pattern = "(PIP2)|(PIP3)", x = swissList$PIP$`Abbreviation*`, ignore.case = FALSE, invert = TRUE)
       # swissList$PIP        = swissList$PIP[idx, ]
       idx                  = grep(pattern = "^PIP-$", x = allAbbrevs, ignore.case = FALSE)
       swissList$PIP        = swissdb[idx, ]
       
       # filter for PIP2
       idx                  = grep(pattern = "^PIP2-$", x = allAbbrevs, ignore.case = FALSE)
       swissList$PIP2       = swissdb[idx, ]
       
       # filter for PIP3
       idx                  = grep(pattern = "^PIP3-$", x = allAbbrevs, ignore.case = FALSE)
       swissList$PIP3       = swissdb[idx, ]
       
       
       # filter for FA
       idx                  = grep(pattern = "^FA-$", x = allAbbrevs, ignore.case = FALSE)
       swissList$FA         = swissdb[idx, ]
       
       
       # filter for CE
       # idx                  = grep(pattern = "(CE)", x = swissdb$`Abbreviation*`, ignore.case = FALSE)
       # swissList$LPS        = swissdb[idx, ]
       
       
       
       
       
       
       return(list(swissList = swissList, 
                   lipidHits = lipidHits, 
                   allClasses = unique(allAbbrevs)))
       
}

