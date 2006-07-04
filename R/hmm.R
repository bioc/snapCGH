#criteria used below are:
	# AIC or BIC
	#if BIC is choosen then it is possible to enter a value for delta as well
##"fitHMM.old" <-
##function (acghList, vr = 0.01, maxiter = 100, criteria = "AIC", delta = NA, full.output = FALSE) 
##{
  #some clunky code so you can put the Criteria argument in characters and still perform the
  #boolean opperators on it below:
##          if( criteria == "AIC") {crit = 1}
##          else if (criteria == "BIC") {crit = 2}
##          else crit = 0
##  
##    if ((crit == 1) || (crit == 2)) {
##    datainfo = acghList$genes
##    dat = log2ratios(acghList)
##    chrom.uniq <- unique(datainfo$Chrom)
##    nstates <- matrix(NA, nrow = length(chrom.uniq), ncol = ncol(dat))

    #matrix template.  It saves me having to define a new empty matrix 6 times in next bit of code
##    	template = matrix(NA,nrow(dat),ncol(dat),dimnames=dimnames(dat))
    #we use template here
##    if (full.output == TRUE) {
# #     segList <- list(M.predicted=template,dispersion=template,state=template,rpred=template,prob=template)
##    }
##    else {
##      segList <- list(M.predicted=template,dispersion=template,state=template,)
##    }
##    if (criteria == "BIC") {
##        if (is.na(delta)) {
##            delta <- c(1)
##        }
##		}
##    for (i in 1:ncol(dat)) {
##        cat("sample is ", i, "  Chromosomes: ")
##		counter = 0  #counter to mark place in segList so we know where to put the values for the next chromosome
##        for (j in 1:length(chrom.uniq)) {
##            cat(j, " ")
##            res <- try(states.hmm.func.old(sample = i, chrom = j, 
##                dat = dat, datainfo = datainfo, vr = vr, maxiter = maxiter, 
##                criteria = crit, delta = delta))
##                  nstates[j, i] <- res$nstates.list
## 				foo = dat[datainfo$Chrom == j,i]	
##				segList$M.predicted[(counter+1):(counter+length(foo)),i] = as.matrix(res$out.list[,1])
## 				segList$dispersion[(counter+1):(counter+length(foo)),i] = as.matrix(res$out.list[,2])
## 		#		segList$obs[(counter+1):(counter+length(foo)),i] = as.matrix(res$out.list[,6])
##                                segList$state[(counter+1):(counter+length(foo)),i] = as.matrix(res$out.list[,3])
# #                               if (full.output == TRUE) { #adding the additional output
##                                  
##                                  segList$rpred[(counter+1):(counter+length(foo)),i] = as.matrix(res$out.list[,4])
##                                  segList$prob[(counter+1):(counter+length(foo)),i] = as.matrix(res$out.list[,5])
##                                }
##           		counter = counter + length(foo)
##        }
##        cat("\n")
##    }
##                segList$M.observed = acghList$M.observed
##                segList$num.states = nstates
##                colnames(segList$num.states) <- colnames(dat)
##		segList$genes <- datainfo[,2:3]
##		new("SegList",segList)
##  }
##        else {
##          cat("You must enter AIC or BIC for the criteria argument\n")
##        }
##          
##}
  
##"states.hmm.func.old" <-
##function (sample, chrom, dat, datainfo = genes, vr = 0.01, 
##    maxiter = 100, criteria = 1, delta = 1) 
##{
##    obs <- dat[datainfo$Chrom == chrom, sample]
##    kb <- datainfo$kb[datainfo$Chrom == chrom]
##    obs.ord <- obs[order(kb)]
##    kb.ord <- kb[order(kb)]
##    ind.nonna <- which(!is.na(obs.ord))
##    y <- obs.ord[ind.nonna]
##    kb <- kb.ord[ind.nonna]
##    numobs <- length(y)
##    if (numobs > 2) 
##        pam2 <- kmeans(y, 2)
##    else pam2 <- 1
##    if (numobs > 3) 
##        pam3 <- kmeans(y, 3)
##    else pam3 <- 1
##    if (numobs > 4) 
##        pam4 <- kmeans(y, 4)
##    else pam4 <- 1
##    if (numobs > 5) 
##        pam5 <- kmeans(y, 5)
##   else pam5 <- 1
##    if (numobs > 2) 
##        mu2 <- c(pam2$centers)
##    else mu2 <- y
##    if (numobs > 3) 
##        mu3 <- c(pam3$centers)
##    else mu3 <- y
##    if (numobs > 4) 
##        mu4 <- c(pam4$centers)
##    else mu4 <- y
##    if (numobs > 5) 
##        mu5 <- c(pam5$centers)
##    else mu5 <- y
##    gamma2 <- matrix(c(0.9, 0.1, 0.1, 0.9), ncol = 2, b = TRUE)
##    gamma3 <- matrix(c(0.9, 0.05, 0.05, 0.05, 0.9, 0.05, 0.05, 
##        0.05, 0.9), ncol = 3, b = TRUE)
##    gamma4 <- matrix(c(0.9, rep(0.1/3, 3), 0.1/3, 0.9, rep(0.1/3, 
##        2), rep(0.1/3, 2), 0.9, 0.1/3, rep(0.1/3, 3), 0.9), ncol = 4, 
##        b = TRUE)
##    gamma5 <- matrix(c(0.9, rep(0.025, 4), 0.025, 0.9, rep(0.025, 
##        3), rep(0.025, 2), 0.9, rep(0.025, 2), rep(0.025, 3), 
##        0.9, 0.025, rep(0.025, 4), 0.9), ncol = 5, b = TRUE)
##    k1 <- 2
##    k2 <- 5
##    k3 <- 10
##    k4 <- 17
##    k5 <- 26
##    z1 <- -sum(log(dnorm(y, mean = mean(y), sd = sqrt(var(y)))))
##    z2 <- try(hidden(y, dist = "normal", cmu = mu1.func, pcmu = mu2, 
##        pshape = vr, pgamma = gamma2, iterlim = maxiter))
##    z3 <- try(hidden(y, dist = "normal", cmu = mu1.func, pcmu = mu3, 
##        pshape = vr, pgamma = gamma3, iterlim = maxiter))
##    z4 <- try(hidden(y, dist = "normal", cmu = mu1.func, pcmu = mu4, 
##        pshape = vr, pgamma = gamma4, iterlim = maxiter))
##    z5 <- try(hidden(y, dist = "normal", cmu = mu1.func, pcmu = mu5, 
##        pshape = vr, pgamma = gamma5, iterlim = maxiter))
##   options(show.error.messages = TRUE)
   
##    if (length(names(z2)) == 0) {
##        z2$maxlik <- NA
##    }
##    if (length(names(z3)) == 0) {
##        z3$maxlik <- NA
##    }
##    if (length(names(z4)) == 0) {
##        z4$maxlik <- NA
##    }
##    if (length(names(z5)) == 0) {
##        z5$maxlik <- NA
##   }
##		if (criteria == 1) {
##            factor <- 2
##        }
##        else if (criteria == 2) {
##                factor <- log(numobs) * delta
##        }
##        lik <- c((z1 + k1 * factor/2), (z2$maxlik + k2 * factor/2), 
##            (z3$maxlik + k3 * factor/2), (z4$maxlik + k4 * factor/2), 
##            (z5$maxlik + k5 * factor/2))
##        likmin <- which.min(lik)

##		switch(likmin, z <- z1, z <- z2, z <- z3, z <- z4, z <- z5)
##		nstates <- likmin

##        if (nstates > 1) {
##            maxstate <- apply(z$filter, 2, which.max)
##            rpred <- z$rpred
##            prob <- apply(z$filter, 2, max)
##            maxstate.unique <- unique(maxstate)
##            pred <- rep(0, length(y))
##           disp <- rep(0, length(y))
##            for (m in 1:length(maxstate.unique)) {
##                pred[maxstate == maxstate.unique[m]] <- median(y[maxstate == 
##                  maxstate.unique[m]])
##                disp[maxstate == maxstate.unique[m]] <- mad(y[maxstate == 
##                  maxstate.unique[m]])
##            }
##        }
##        else {
##            maxstate <- rep(1, length(y))
##            rpred <- rep(median(y), length(y))
##            prob <- rep(1, length(y))
##            pred <- rep(median(y), length(y))
##            disp <- rep(mad(y), length(y))
##        }
##        out <- cbind(matrix(pred, ncol = 1), matrix(disp, ncol = 1),matrix(maxstate, ncol = 1), matrix(rpred, 
##            ncol = 1), matrix(prob, ncol = 1))
##        out.all <- matrix(NA, nrow = length(kb.ord), ncol = 6)
##       out.all[ind.nonna, 1:5] <- out
##        out.all[, 6] <- obs.ord
##        out.all <- as.data.frame(out.all)
##        dimnames(out.all)[[2]] <- c("pred", "disp", "state", 
##            "rpred", "prob", "obs")
##    list(out.list = out.all, nstates.list = nstates)
##}

"runHomHMM" <- 
function (MA, vr = 0.01, maxiter = 100, criteria = "AIC", delta = NA, full.output = FALSE, eps = 0.01) 
{

  if (is.null(MA$design)) 
        stop("MA$design component is null")

  for(i in 1:length(MA$design)){
  temp <- MA$design[i]* MA$M[,i]
  MA$M[,i] <- temp
  }
##Changing the names of the Chr and Position columns so that the aCGH code can access them.

  colnames(MA$genes)[colnames(MA$genes) == "Position"] = "kb"
  colnames(MA$genes)[colnames(MA$genes) == "Chr"] = "Chrom"

  #colnames(MA$genes)[[6]] <- "kb"
  #colnames(MA$genes)[[7]] <- "Chrom"

  #some clunky code so you can put the Criteria argument in characters and still perform the
  #boolean opperators on it below:
  crit = TRUE
  if( criteria == "AIC") {
	aic = TRUE
	bic = FALSE
	}
  else if (criteria == "BIC") {
	bic = TRUE
	aic = FALSE
	}
  else crit = FALSE
  
    if (crit == TRUE) {
    datainfo = MA$genes
    dat = log2ratios(MA)
    chrom.uniq <- unique(datainfo$Chrom)
    nstates <- matrix(NA, nrow = length(chrom.uniq), ncol = ncol(dat))

    #matrix template.  It saves me having to define a new empty matrix 6 times in next bit of code
    	template = matrix(NA,nrow(dat),ncol(dat),dimnames=dimnames(dat))


    #we use template here
    if (full.output == TRUE) {
      segList <- list(M.predicted=template,dispersion=template,state=template,smoothed=template,probability=template)
    }
    else {
      segList <- list(M.predicted=template,state=template)
    }
    if (criteria == "BIC") {
        if (is.na(delta)) {
            delta <- c(1)
        }
		}
    for (i in 1:ncol(dat)) {
        cat("sample is ", i, "  Chromosomes: ")
		counter = 0  #counter to mark place in segList so we know where to put the values for the next chromosome
        for (j in 1:length(chrom.uniq)) {
            cat(chrom.uniq[j], " ")
            foo = dat[datainfo$Chrom == chrom.uniq[j],i]
            if(length(foo) > 5){
              res <- try(states.hmm.func(sample = i, chrom = chrom.uniq[j], 
                dat = dat, datainfo = datainfo, vr = vr, maxiter = maxiter, 
                aic = aic, bic = bic, delta = delta, eps = eps, nlists = 1))
                  nstates[j, i] <- res$nstates.list[[1]]
 				segList$M.predicted[(counter+1):(counter+length(foo)),i] = as.matrix(res$out.list[[1]]$pred)			
                                segList$state[(counter+1):(counter+length(foo)),i] = as.matrix(res$out.list[[1]]$state)
                                if (full.output == TRUE) { #adding the additional output.
                                  segList$dispersion[(counter+1):(counter+length(foo)),i] = as.matrix(res$out.list[[1]]$disp)
                                  segList$smoothed[(counter+1):(counter+length(foo)),i] = as.matrix(res$out.list[[1]]$rpred)
                                  segList$probability[(counter+1):(counter+length(foo)),i] = as.matrix(res$out.list[[1]]$prob)
                                }
            		counter = counter + length(foo)
          }
            else {
              cat("\nToo few observations.  See help file for how this is handled\n")
              nstates[j, i] <- 1
              segList$M.predicted[(counter+1):(counter+length(foo)),i] = rep(mean(foo), length(foo))
              segList$state[(counter+1):(counter+length(foo)),i] = rep(1, length(foo))
              counter = counter + length(foo)
            }
        }
        cat("\n")
    }
                segList$M.observed = MA$M
                segList$num.states = nstates
                colnames(segList$num.states) <- colnames(dat)
                rownames(segList$num.states) <- paste("Chrom",unique(MA$genes$Chr))
		segList$genes <- datainfo
    
                ##Changing the names of the Chr and Position back.

    colnames(segList$genes)[colnames(segList$genes) == "kb"] = "Position"
    colnames(segList$genes)[colnames(segList$genes) == "Chrom"] = "Chr"
#                colnames(segList$genes)[[6]] <- "Position"
#                colnames(segList$genes)[[7]] <- "Chr"

		new("SegList",segList)
  }
        else {
          cat("You must enter AIC or BIC for the criteria argument\n")
        }
          
}

