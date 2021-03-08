# for (lipidSpecies in names(searchList$lipidHits))
# {
#        
#        
#        if(nrow(searchList$lipidHits[[lipidSpecies]]) == 0) {next}
#        
#        par(mfrow = c(2,2))
#        
#        hitsppp               = spatstat::ppp(x = searchList$lipidHits[[lipidSpecies]]$x,
#                                              y = searchList$lipidHits[[lipidSpecies]]$y,
#                                              window = winList[[regName]],
#                                              marks = searchList$lipidHits[[lipidSpecies]][ , -c(1, 2)])
#        
#        hitsK                = spatstat::subset.ppp(hitsppp, modeAdduct == "M+K")
#        hitsNa               = spatstat::subset.ppp(hitsppp, modeAdduct == "M+Na")
#        
#        if(hitsK$n == 0 || hitsNa$n == 0) {next}
#        
#        #hitsppp$marks$intensity = sqrt(hitsppp$marks$intensity)
#        # hitsppp$marks$intensity = (hitsppp$marks$intensity - min(hitsppp$marks$intensity)) /
#        #        (max(hitsppp$marks$intensity) - min(hitsppp$marks$intensity))
#        
#        
#        #// H&E ----
#        # EBImage::display(img1, method = "raster")
#        # legend("topleft", legend = c("MVP", "CT", "IT", "LE", "Necrosis"), 
#        #               col = c("red", "blue", "green", "yellow", "black"),
#        #               bty = "n", horiz = TRUE, pch = 19)
#        #    
#        
#        
#        #//_____________________________________________________________________________ K+ ----
#        pppPixalateAllK      = spatstat::pixellate(hitsK, 
#                                                   weights = hitsK$marks$intensity, 
#                                                   W = spatstat::as.mask(winList[[regName]], 
#                                                                         dimyx=c(diff(winList[[regName]]$yrange) + 1, 
#                                                                                 diff(winList[[regName]]$xrange) + 1)), 
#                                                   padzero = FALSE, savemap = TRUE)
#        
#        spatstat::plot.im(pppPixalateAllK, 
#                          main = paste0(lipidSpecies," - ",length(unique(hitsK$marks$mass)), 
#                                        " K+ adduct formations detected"), 
#                          ylim = rev(range(winList[[regName]]$y)))
#        
#        # annotations
#        spatstat::plot.owin(annotationList$MVP, border = "red", lty = "dashed", lwd = 1.5, add = TRUE)
#        spatstat::plot.owin(annotationList$CT, border = "blue",lty = "dashed" , lwd = 1.5, add = TRUE)
#        spatstat::plot.owin(annotationList$IT, border = "green",lty = "dashed" , lwd = 1.5, add = TRUE)
#        spatstat::plot.owin(annotationList$LE, border = "yellow",lty = "dashed" , lwd = 1.5, add = TRUE)
#        spatstat::plot.owin(annotationList$nec, border = "black",lty = "dashed" , lwd = 1.5, add = TRUE)
#        
#        
#        #//_____________________________________________________________________________ K+ significant ----
#        
#        ## craete a complete spatial randomness point pattern with the same number of points and window
#        csrppp               = spatstat::rpoispp(lambda = spatstat::intensity(hitsK), win = winList[[regName]]) # sigma = 3.2245
#        csrppp$marks         = data.frame(intensity = rpois(length(csrppp$n), mean(hitsK$marks$intensity)))
#        #spatstat::intensity(subsetppp)
#        
#        # create a density map for csr
#        #sppp                 = spatstat::bw.ppl(unique(csrppp, rule="unmark", warn = FALSE), shortcut=TRUE) # computes sigma bases on spp
#        sppp                 = 3.2245
#        
#        ##unique(subsetppp, rule="unmark", warn = FALSE)
#        denCsrAllK           = spatstat::density.ppp(x = csrppp, sigma = sppp, weights = csrppp$marks$intensity)
#        
#        # calculate the probability function Fxy by normalizing denHeme to denCsr
#        Fxy                  = (pppPixalateAllK - mean(denCsrAllK)) / sd(denCsrAllK)
#        
#        # cut-off based on inverse standrard normal cumulative density function
#        cutoff               = qnorm(0.05, lower.tail = F)
#        
#        # with Bonferroni correction
#        cutoff               = qnorm(0.05/hitsK$n, lower.tail = F)
#        
#        
#        # normalized probability function - Kather et al. 
#        # hotspotIm            = spatstat::eval.im(Fxy * (Fxy >= cutoff))
#        # hotspotIm[hotspotIm == 0] = NA
#        # spatstat::plot.im(hotspotIm, clipwin = winList[[regName]], 
#        #             main = paste0("Normalized prob. fun. - Phosphatidylserine "), 
#        #             ylim = rev(winList[[regName]]$yrange), frame.plot = FALSE)
#        # spatstat::plot.owin(winList[[regName]], add = TRUE)      
#        
#        
#        # combine the two weighted figures 
#        spatstat::plot.im(pppPixalateAllK, 
#                          main = paste0(lipidSpecies,"-K+ significant increase overlayed-all"), 
#                          ylim = rev(range(winList[[regName]]$y)))
#        
#        
#        
#        hotspotIm            = spatstat::eval.im(Fxy * (Fxy >= cutoff))
#        hotspotIm[hotspotIm > 0] = NA
#        
#        spatstat::plot.im(hotspotIm, col = rgb(0,0,0,0.7), #spatstat::colourmap(col = rgb(0,0,0, 0.5), range = range(hotspotIm))
#                          clipwin = winList[[regName]], 
#                          frame.plot = FALSE, add = TRUE)
#        
#        # annotations
#        spatstat::plot.owin(annotationList$MVP, border = "red", lty = "dashed", lwd = 1.5, add = TRUE)
#        spatstat::plot.owin(annotationList$CT, border = "blue",lty = "dashed" , lwd = 1.5, add = TRUE)
#        spatstat::plot.owin(annotationList$IT, border = "green",lty = "dashed" , lwd = 1.5, add = TRUE)
#        spatstat::plot.owin(annotationList$LE, border = "yellow",lty = "dashed" , lwd = 1.5, add = TRUE)
#        spatstat::plot.owin(annotationList$nec, border = "black",lty = "dashed" , lwd = 1.5, add = TRUE)
#        
#        #//_____________________________________________________________________________ K+ confirmed ----
#        
#        hitsKConf            = spatstat::subset.ppp(hitsK, mass %in% confirmedHistList[[lipidSpecies]])
#        
#        pppPixalateConfK      = spatstat::pixellate(hitsKConf, 
#                                                    weights = hitsKConf$marks$intensity, 
#                                                    W = spatstat::as.mask(winList[[regName]], 
#                                                                          dimyx=c(diff(winList[[regName]]$yrange) + 1, 
#                                                                                  diff(winList[[regName]]$xrange) + 1)), 
#                                                    padzero = FALSE, savemap = TRUE)
#        
#        spatstat::plot.im(pppPixalateConfK, 
#                          main = paste0(lipidSpecies," - ",length(unique(hitsKConf$marks$mass)), 
#                                        " unique K+ adduct formations confirmed"), 
#                          ylim = rev(range(winList[[regName]]$y)))
#        
#        # annotations
#        spatstat::plot.owin(annotationList$MVP, border = "red", lty = "dashed", lwd = 1.5, add = TRUE)
#        spatstat::plot.owin(annotationList$CT, border = "blue",lty = "dashed" , lwd = 1.5, add = TRUE)
#        spatstat::plot.owin(annotationList$IT, border = "green",lty = "dashed" , lwd = 1.5, add = TRUE)
#        spatstat::plot.owin(annotationList$LE, border = "yellow",lty = "dashed" , lwd = 1.5, add = TRUE)
#        spatstat::plot.owin(annotationList$nec, border = "black",lty = "dashed" , lwd = 1.5, add = TRUE)
#        
#        
#        #//_____________________________________________________________________________ K+ confirmed significant ----
#        
#        ## craete a complete spatial randomness point pattern with the same number of points and window
#        csrppp               = spatstat::rpoispp(lambda = spatstat::intensity(hitsKConf), win = winList[[regName]]) # sigma = 3.2245
#        csrppp$marks         = data.frame(intensity = rpois(length(csrppp$n), mean(hitsKConf$marks$intensity)))
#        #spatstat::intensity(subsetppp)
#        
#        # create a density map for csr
#        #sppp                 = spatstat::bw.ppl(unique(csrppp, rule="unmark", warn = FALSE), shortcut=TRUE) # computes sigma bases on spp
#        sppp                 = 3.2245
#        
#        ##unique(subsetppp, rule="unmark", warn = FALSE)
#        denCsrConfK          = spatstat::density.ppp(x = csrppp, sigma = sppp, weights = csrppp$marks$intensity)
#        
#        # calculate the probability function Fxy by normalizing denHeme to denCsr
#        Fxy                  = (pppPixalateConfK - mean(denCsrConfK)) / sd(denCsrConfK)
#        
#        # cut-off based on inverse standrard normal cumulative density function
#        cutoff               = qnorm(0.05, lower.tail = F)
#        
#        # with Bonferroni correction
#        cutoff               = qnorm(0.05/hitsKConf$n, lower.tail = F)
#        
#        
#        # normalized probability function - Kather et al. 
#        # hotspotIm            = spatstat::eval.im(Fxy * (Fxy >= cutoff))
#        # hotspotIm[hotspotIm == 0] = NA
#        # spatstat::plot.im(hotspotIm, clipwin = winList[[regName]], 
#        #             main = paste0("Normalized prob. fun. - Phosphatidylserine "), 
#        #             ylim = rev(winList[[regName]]$yrange), frame.plot = FALSE)
#        # spatstat::plot.owin(winList[[regName]], add = TRUE)      
#        
#        
#        # combine the two weighted figures 
#        spatstat::plot.im(pppPixalateConfK, 
#                          main = paste0(lipidSpecies,"-K+ significant increase overlayed-confirmed"), 
#                          ylim = rev(range(winList[[regName]]$y)))
#        
#        
#        
#        hotspotIm            = spatstat::eval.im(Fxy * (Fxy >= cutoff))
#        hotspotIm[hotspotIm > 0] = NA
#        
#        spatstat::plot.im(hotspotIm, col = rgb(0,0,0,0.7), #spatstat::colourmap(col = rgb(0,0,0, 0.5), range = range(hotspotIm))
#                          clipwin = winList[[regName]], 
#                          frame.plot = FALSE, add = TRUE)
#        
#        # annotations
#        spatstat::plot.owin(annotationList$MVP, border = "red", lty = "dashed", lwd = 1.5, add = TRUE)
#        spatstat::plot.owin(annotationList$CT, border = "blue",lty = "dashed" , lwd = 1.5, add = TRUE)
#        spatstat::plot.owin(annotationList$IT, border = "green",lty = "dashed" , lwd = 1.5, add = TRUE)
#        spatstat::plot.owin(annotationList$LE, border = "yellow",lty = "dashed" , lwd = 1.5, add = TRUE)
#        spatstat::plot.owin(annotationList$nec, border = "black",lty = "dashed" , lwd = 1.5, add = TRUE)
#        
#        
#        
#        
#        
#        #//_____________________________________________________________________________ Na+ ----
#        pppPixalateAllNa      = spatstat::pixellate(hitsNa, 
#                                                    weights = hitsNa$marks$intensity, 
#                                                    W = spatstat::as.mask(winList[[regName]], 
#                                                                          dimyx=c(diff(winList[[regName]]$yrange) + 1, 
#                                                                                  diff(winList[[regName]]$xrange) + 1)), 
#                                                    padzero = FALSE, savemap = TRUE)
#        
#        spatstat::plot.im(pppPixalateAllNa, 
#                          main = paste0(lipidSpecies," - ",length(unique(hitsNa$marks$mass)), 
#                                        " Na+ adduct formations detected"), 
#                          ylim = rev(range(winList[[regName]]$y)))
#        
#        # annotations
#        spatstat::plot.owin(annotationList$MVP, border = "red", lty = "dashed", lwd = 1.5, add = TRUE)
#        spatstat::plot.owin(annotationList$CT, border = "blue",lty = "dashed" , lwd = 1.5, add = TRUE)
#        spatstat::plot.owin(annotationList$IT, border = "green",lty = "dashed" , lwd = 1.5, add = TRUE)
#        spatstat::plot.owin(annotationList$LE, border = "yellow",lty = "dashed" , lwd = 1.5, add = TRUE)
#        spatstat::plot.owin(annotationList$nec, border = "black",lty = "dashed" , lwd = 1.5, add = TRUE)
#        
#        
#        #//_____________________________________________________________________________ Na+ significant ----
#        
#        ## craete a complete spatial randomness point pattern with the same number of points and window
#        csrppp               = spatstat::rpoispp(lambda = spatstat::intensity(hitsNa), win = winList[[regName]]) # sigma = 3.2245
#        csrppp$marks         = data.frame(intensity = rpois(length(csrppp$n), mean(hitsNa$marks$intensity)))
#        #spatstat::intensity(subsetppp)
#        
#        # create a density map for csr
#        #sppp                 = spatstat::bw.ppl(unique(csrppp, rule="unmark", warn = FALSE), shortcut=TRUE) # computes sigma bases on spp
#        sppp                 = 3.2245
#        
#        ##unique(subsetppp, rule="unmark", warn = FALSE)
#        denCsrAllNa          = spatstat::density.ppp(x = csrppp, sigma = sppp, weights = csrppp$marks$intensity)
#        
#        # calculate the probability function Fxy by normalizing denHeme to denCsr
#        Fxy                  = (pppPixalateAllNa - mean(denCsrAllNa)) / sd(denCsrAllNa)
#        
#        # cut-off based on inverse standrard normal cumulative density function
#        cutoff               = qnorm(0.05, lower.tail = F)
#        
#        # with Bonferroni correction
#        cutoff               = qnorm(0.05/hitsNa$n, lower.tail = F)
#        
#        
#        # normalized probability function - Kather et al. 
#        # hotspotIm            = spatstat::eval.im(Fxy * (Fxy >= cutoff))
#        # hotspotIm[hotspotIm == 0] = NA
#        # spatstat::plot.im(hotspotIm, clipwin = winList[[regName]], 
#        #             main = paste0("Normalized prob. fun. - Phosphatidylserine "), 
#        #             ylim = rev(winList[[regName]]$yrange), frame.plot = FALSE)
#        # spatstat::plot.owin(winList[[regName]], add = TRUE)      
#        
#        
#        # combine the two weighted figures 
#        spatstat::plot.im(pppPixalateAllNa, 
#                          main = paste0(lipidSpecies,"-Na+ significant increase overlayed-all"), 
#                          ylim = rev(range(winList[[regName]]$y)))
#        
#        
#        
#        hotspotIm            = spatstat::eval.im(Fxy * (Fxy >= cutoff))
#        hotspotIm[hotspotIm > 0] = NA
#        
#        spatstat::plot.im(hotspotIm, col = rgb(0,0,0,0.7), #spatstat::colourmap(col = rgb(0,0,0, 0.5), range = range(hotspotIm))
#                          clipwin = winList[[regName]], 
#                          frame.plot = FALSE, add = TRUE)
#        
#        # annotations
#        spatstat::plot.owin(annotationList$MVP, border = "red", lty = "dashed", lwd = 1.5, add = TRUE)
#        spatstat::plot.owin(annotationList$CT, border = "blue",lty = "dashed" , lwd = 1.5, add = TRUE)
#        spatstat::plot.owin(annotationList$IT, border = "green",lty = "dashed" , lwd = 1.5, add = TRUE)
#        spatstat::plot.owin(annotationList$LE, border = "yellow",lty = "dashed" , lwd = 1.5, add = TRUE)
#        spatstat::plot.owin(annotationList$nec, border = "black",lty = "dashed" , lwd = 1.5, add = TRUE)
#        
#        #//_____________________________________________________________________________ Na+ confirmed ----
#        
#        hitsNaConf            = spatstat::subset.ppp(hitsNa, mass %in% confirmedHistList[[lipidSpecies]])
#        
#        pppPixalateConfNa      = spatstat::pixellate(hitsNaConf, 
#                                                     weights = hitsNaConf$marks$intensity, 
#                                                     W = spatstat::as.mask(winList[[regName]], 
#                                                                           dimyx=c(diff(winList[[regName]]$yrange) + 1, 
#                                                                                   diff(winList[[regName]]$xrange) + 1)), 
#                                                     padzero = FALSE, savemap = TRUE)
#        
#        spatstat::plot.im(pppPixalateConfNa, 
#                          main = paste0(lipidSpecies," - ",length(unique(hitsNaConf$marks$mass)), 
#                                        " unique Na+ adduct formations confirmed"), 
#                          ylim = rev(range(winList[[regName]]$y)))
#        
#        # annotations
#        spatstat::plot.owin(annotationList$MVP, border = "red", lty = "dashed", lwd = 1.5, add = TRUE)
#        spatstat::plot.owin(annotationList$CT, border = "blue",lty = "dashed" , lwd = 1.5, add = TRUE)
#        spatstat::plot.owin(annotationList$IT, border = "green",lty = "dashed" , lwd = 1.5, add = TRUE)
#        spatstat::plot.owin(annotationList$LE, border = "yellow",lty = "dashed" , lwd = 1.5, add = TRUE)
#        spatstat::plot.owin(annotationList$nec, border = "black",lty = "dashed" , lwd = 1.5, add = TRUE)
#        
#        
#        #//_____________________________________________________________________________ Na+ confirmed significant ----
#        
#        ## craete a complete spatial randomness point pattern with the same number of points and window
#        csrppp               = spatstat::rpoispp(lambda = spatstat::intensity(hitsNaConf), win = winList[[regName]]) # sigma = 3.2245
#        csrppp$marks         = data.frame(intensity = rpois(length(csrppp$n), mean(hitsNaConf$marks$intensity)))
#        #spatstat::intensity(subsetppp)
#        
#        # create a density map for csr
#        #sppp                 = spatstat::bw.ppl(unique(csrppp, rule="unmark", warn = FALSE), shortcut=TRUE) # computes sigma bases on spp
#        sppp                 = 3.2245
#        
#        ##unique(subsetppp, rule="unmark", warn = FALSE)
#        denCsrConfNa          = spatstat::density.ppp(x = csrppp, sigma = sppp, weights = csrppp$marks$intensity)
#        
#        # calculate the probability function Fxy by normalizing denHeme to denCsr
#        Fxy                  = (pppPixalateConfNa - mean(denCsrConfNa)) / sd(denCsrConfNa)
#        
#        # cut-off based on inverse standrard normal cumulative density function
#        cutoff               = qnorm(0.05, lower.tail = F)
#        
#        # with Bonferroni correction
#        cutoff               = qnorm(0.05/hitsNaConf$n, lower.tail = F)
#        
#        
#        # normalized probability function - Kather et al. 
#        # hotspotIm            = spatstat::eval.im(Fxy * (Fxy >= cutoff))
#        # hotspotIm[hotspotIm == 0] = NA
#        # spatstat::plot.im(hotspotIm, clipwin = winList[[regName]], 
#        #             main = paste0("Normalized prob. fun. - Phosphatidylserine "), 
#        #             ylim = rev(winList[[regName]]$yrange), frame.plot = FALSE)
#        # spatstat::plot.owin(winList[[regName]], add = TRUE)      
#        
#        
#        # combine the two weighted figures 
#        spatstat::plot.im(pppPixalateConfNa, 
#                          main = paste0(lipidSpecies,"-Na+ significant increase overlayed-confirmed"), 
#                          ylim = rev(range(winList[[regName]]$y)))
#        
#        
#        
#        hotspotIm            = spatstat::eval.im(Fxy * (Fxy >= cutoff))
#        hotspotIm[hotspotIm > 0] = NA
#        
#        spatstat::plot.im(hotspotIm, col = rgb(0,0,0,0.7), #spatstat::colourmap(col = rgb(0,0,0, 0.5), range = range(hotspotIm))
#                          clipwin = winList[[regName]], 
#                          frame.plot = FALSE, add = TRUE)
#        
#        # annotations
#        spatstat::plot.owin(annotationList$MVP, border = "red", lty = "dashed", lwd = 1.5, add = TRUE)
#        spatstat::plot.owin(annotationList$CT, border = "blue",lty = "dashed" , lwd = 1.5, add = TRUE)
#        spatstat::plot.owin(annotationList$IT, border = "green",lty = "dashed" , lwd = 1.5, add = TRUE)
#        spatstat::plot.owin(annotationList$LE, border = "yellow",lty = "dashed" , lwd = 1.5, add = TRUE)
#        spatstat::plot.owin(annotationList$nec, border = "black",lty = "dashed" , lwd = 1.5, add = TRUE)
#        
#        
#        
#        
#        #//_____________________________________________________________________________ K+/Na+ all  ----
#        
#        #// K+/Na+
#        ionRatio             = spatstat::eval.im(pppPixalateAllK / pppPixalateAllNa)
#        #upperlim             = quantile(ionRatio$v, 0.95, na.rm = T)
#        #ionRatio$v[ionRatio$v > upperlim] = upperlim
#        
#        spatstat::plot.im(ionRatio, main = paste0(lipidSpecies, "- Density map of K+/Na+"), 
#                          ylim = rev(ionRatio$yrange), log = F)
#        
#        # annotations
#        spatstat::plot.owin(annotationList$MVP, border = "red", lty = "dashed", lwd = 1.5, add = TRUE)
#        spatstat::plot.owin(annotationList$CT, border = "blue",lty = "dashed" , lwd = 1.5, add = TRUE)
#        spatstat::plot.owin(annotationList$IT, border = "green",lty = "dashed" , lwd = 1.5, add = TRUE)
#        spatstat::plot.owin(annotationList$LE, border = "yellow",lty = "dashed" , lwd = 1.5, add = TRUE)
#        spatstat::plot.owin(annotationList$nec, border = "black",lty = "dashed" , lwd = 1.5, add = TRUE)
#        
#        
#        
#        
#        #//_____________________________________________________________________________ K+/Na+ all significant ----
#        
#        
#        ##unique(subsetppp, rule="unmark", warn = FALSE)
#        denCsrRatio          = spatstat::eval.im(denCsrAllK / denCsrAllNa)
#        
#        # calculate the probability function Fxy by normalizing denHeme to denCsr
#        Fxy                  = (ionRatio - mean(denCsrRatio)) / sd(denCsrRatio)
#        
#        # cut-off based on inverse standrard normal cumulative density function
#        cutoff               = qnorm(0.05, lower.tail = F)
#        
#        # with Bonferroni correction
#        cutoff               = qnorm(0.05/hitsK$n, lower.tail = F)
#        
#        
#        # normalized probability function - Kather et al. 
#        # hotspotIm            = spatstat::eval.im(Fxy * (Fxy >= cutoff))
#        # hotspotIm[hotspotIm == 0] = NA
#        # spatstat::plot.im(hotspotIm, clipwin = winList[[regName]], 
#        #             main = paste0("Normalized prob. fun. - Phosphatidylserine "), 
#        #             ylim = rev(winList[[regName]]$yrange), frame.plot = FALSE)
#        # spatstat::plot.owin(winList[[regName]], add = TRUE)      
#        
#        
#        # combine the two weighted figures 
#        spatstat::plot.im(ionRatio, 
#                          main = paste0(lipidSpecies," K+/Na+ significant increase overlayed-all"), 
#                          ylim = rev(range(winList[[regName]]$y)))
#        
#        
#        
#        hotspotIm            = spatstat::eval.im(Fxy * (Fxy >= cutoff))
#        hotspotIm[hotspotIm > 0] = NA
#        
#        spatstat::plot.im(hotspotIm, col = rgb(0,0,0,0.7), #spatstat::colourmap(col = rgb(0,0,0, 0.5), range = range(hotspotIm))
#                          clipwin = winList[[regName]], 
#                          frame.plot = FALSE, add = TRUE)
#        
#        # annotations
#        spatstat::plot.owin(annotationList$MVP, border = "red", lty = "dashed", lwd = 1.5, add = TRUE)
#        spatstat::plot.owin(annotationList$CT, border = "blue",lty = "dashed" , lwd = 1.5, add = TRUE)
#        spatstat::plot.owin(annotationList$IT, border = "green",lty = "dashed" , lwd = 1.5, add = TRUE)
#        spatstat::plot.owin(annotationList$LE, border = "yellow",lty = "dashed" , lwd = 1.5, add = TRUE)
#        spatstat::plot.owin(annotationList$nec, border = "black",lty = "dashed" , lwd = 1.5, add = TRUE)
#        
#        
#        
#        
#        
#        #//_____________________________________________________________________________ K+/Na+ conf  ----
#        
#        #// K+/Na+
#        ionRatio             = spatstat::eval.im(pppPixalateConfK / pppPixalateConfNa)
#        #upperlim             = quantile(ionRatio$v, 0.95, na.rm = T)
#        #ionRatio$v[ionRatio$v > upperlim] = upperlim
#        
#        
#        
#        spatstat::plot.im(ionRatio, main = paste0(lipidSpecies, "- Density map of K+/Na+ confirmed"), 
#                          ylim = rev(ionRatio$yrange), log = F)
#        
#        # annotations
#        spatstat::plot.owin(annotationList$MVP, border = "red", lty = "dashed", lwd = 1.5, add = TRUE)
#        spatstat::plot.owin(annotationList$CT, border = "blue",lty = "dashed" , lwd = 1.5, add = TRUE)
#        spatstat::plot.owin(annotationList$IT, border = "green",lty = "dashed" , lwd = 1.5, add = TRUE)
#        spatstat::plot.owin(annotationList$LE, border = "yellow",lty = "dashed" , lwd = 1.5, add = TRUE)
#        spatstat::plot.owin(annotationList$nec, border = "black",lty = "dashed" , lwd = 1.5, add = TRUE)
#        
#        
#        #//_____________________________________________________________________________ K+/Na+ conf significant ----
#        
#        
#        ##unique(subsetppp, rule="unmark", warn = FALSE)
#        denCsrRatio          = spatstat::eval.im(denCsrConfK / denCsrConfNa)
#        
#        # calculate the probability function Fxy by normalizing denHeme to denCsr
#        Fxy                  = (ionRatio - mean(denCsrRatio)) / sd(denCsrRatio)
#        
#        # cut-off based on inverse standrard normal cumulative density function
#        cutoff               = qnorm(0.05, lower.tail = F)
#        
#        # with Bonferroni correction
#        cutoff               = qnorm(0.05/hitsKConf$n, lower.tail = F)
#        
#        
#        # normalized probability function - Kather et al. 
#        # hotspotIm            = spatstat::eval.im(Fxy * (Fxy >= cutoff))
#        # hotspotIm[hotspotIm == 0] = NA
#        # spatstat::plot.im(hotspotIm, clipwin = winList[[regName]], 
#        #             main = paste0("Normalized prob. fun. - Phosphatidylserine "), 
#        #             ylim = rev(winList[[regName]]$yrange), frame.plot = FALSE)
#        # spatstat::plot.owin(winList[[regName]], add = TRUE)      
#        
#        
#        # combine the two weighted figures 
#        spatstat::plot.im(ionRatio, 
#                          main = paste0(lipidSpecies," K+/Na+ significant increase overlayed-confirmed"), 
#                          ylim = rev(range(winList[[regName]]$y)))
#        
#        
#        
#        hotspotIm            = spatstat::eval.im(Fxy * (Fxy >= cutoff))
#        hotspotIm[hotspotIm > 0] = NA
#        
#        spatstat::plot.im(hotspotIm, col = rgb(0,0,0,0.7), #spatstat::colourmap(col = rgb(0,0,0, 0.5), range = range(hotspotIm))
#                          clipwin = winList[[regName]], 
#                          frame.plot = FALSE, add = TRUE)
#        
#        # annotations
#        spatstat::plot.owin(annotationList$MVP, border = "red", lty = "dashed", lwd = 1.5, add = TRUE)
#        spatstat::plot.owin(annotationList$CT, border = "blue",lty = "dashed" , lwd = 1.5, add = TRUE)
#        spatstat::plot.owin(annotationList$IT, border = "green",lty = "dashed" , lwd = 1.5, add = TRUE)
#        spatstat::plot.owin(annotationList$LE, border = "yellow",lty = "dashed" , lwd = 1.5, add = TRUE)
#        spatstat::plot.owin(annotationList$nec, border = "black",lty = "dashed" , lwd = 1.5, add = TRUE)
#        
#        
#        
#        
#        
# }