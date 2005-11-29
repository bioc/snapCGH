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
##    dat = log2.ratios(acghList)
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



states.hmm.func <- function(sample, chrom, dat, datainfo = clones.info, vr = .01,
             maxiter = 100, criteria = 1, delta = 1,
             nlists = 1, eps = .01, print.info = FALSE,
             diag.prob = .99)
{

    obs <- dat[datainfo$Chr==chrom, sample]
    kb <- datainfo$Position[datainfo$Chr==chrom]
    ##with current sproc files, data is already ordered by kb's
    obs.ord <- obs[order(kb)]
    kb.ord <- kb[order(kb)]

    ind.nonna <- which(!is.na(obs.ord))

    y <- obs.ord[ind.nonna]
    kb <- kb.ord[ind.nonna]


#####################################

    numobs <- length(y)
    zz <- vector(mode = "list", 5)
    zz[[1]] <-
        list(log.lik =
             sum(dnorm(y, mean = mean(y), sd = sd(y), log = T)))
    for(k in 2:5)
    {
##Code added by Mike to make this work if there are less than 5 observations
      if(numobs > k){      
        mu <- kmeans(y, k)$centers
      }
      else{
        mu <- y
      }
#      mu <- kmeans(y, k)$centers
##Finished editing here
      gamma <- matrix((1 - diag.prob) / (k - 1), k, k)
        diag(gamma) <- diag.prob
        zz[[k]] <-
        {

            res <-
                .C("calc_observed_likelihood",
                   as.integer(numobs),
                   as.double(y),
                   as.integer(k),
                   mu = as.double(mu),
                   sigma = as.double(sqrt(vr)),
                   gamma = as.double(gamma),
                   pi = as.double(rep(-log(k), k)),
                   num.iter = as.integer(maxiter),
                   as.double(eps),
                   log.lik = double(1),
                   filtered.cond.probs = double(k * numobs),
                   hidden.states = integer(numobs),
                   as.logical(print.info),
                   PACKAGE = "snapCGH")
            res$hidden.states <- res$hidden.states + 1
            res$filtered.cond.probs <-
                matrix(res$filtered.cond.probs, nr = k)
            res$gamma <- exp(matrix(res$gamma, nr = k))
            res
            
        }
        
    }

###############################################3
###############################################3
    ##identify the model with the smallest model selection criteria

    ##now, scroll over all options:

    ##number of states (means) + number of states*(number of states-1) (transitions) #+ 1 (variance)

#    kk <- c(2, 5, 10, 17, 26)
    kk <- (1:5) ^ 2 + 1
    for (nl in 1:nlists)
    {
      if (criteria == 1)
        {
          factor <- 2
        }
      else if (criteria == 2)
        {
          factor <- log(numobs)*delta
        }
    
        lik <- sapply(zz, function(z) -z$log.lik) + kk * factor / 2
        nstates <- likmin <- which.min(lik)
        z <- zz[[likmin]]

######################################
        ##out rpred and state

        if (nstates > 1) #if non-generic
        {
            ##print(nstates)
            maxstate <- apply(z$filter, 2, which.max)
###            maxstate <- z$hidden.states
            rpred <- as.vector(z$mu %*% z$filter)
            prob <- apply(z$filter, 2, max)
            ##use median for prediction and mad for state dispersions
            maxstate.unique <- unique(maxstate)
            pred <- rep(0, length(y))
            disp <- rep(0, length(y))
            for (m in 1:length(maxstate.unique))
            {
                
                pred[maxstate==maxstate.unique[m]] <-
                    median(y[maxstate==maxstate.unique[m]])
                disp[maxstate==maxstate.unique[m]] <-
                    mad(y[maxstate==maxstate.unique[m]])
                
            }

            ##if (length(z$pshape) == 1)
            ##{
            ##        disp <- rep(z$pshape, length(maxstate))
            ##}
            ##else
            ##{
            ##        disp <- z$pshape[maxstate]
            ##}
            
        }
        else #if generic
        {
            maxstate <- rep(1, length(y))
            ##rpred <- rep(mean(y), length(y))
            rpred <- rep(median(y), length(y))
            prob <- rep(1, length(y))
            ##pred <- rep(mean(y), length(y))
            pred <- rep(median(y), length(y))
            ##disp <- rep(var(y), length(y))
            disp <- rep(mad(y), length(y))
            
        }
        
        out <-
            cbind(matrix(maxstate, ncol=1),
                  matrix(rpred, ncol=1),
                  matrix(prob, ncol=1),
                  matrix(pred, ncol=1),
                  matrix(disp, ncol=1))
        
        out.all <- matrix(NA, nrow=length(kb.ord), ncol=6)
        out.all[ind.nonna,1:5] <- out
        
        out.all[,6] <- obs.ord
        out.all <- as.data.frame(out.all)
        dimnames(out.all)[[2]] <- c("state", "rpred", "prob", "pred", "dispersion", "M.observed")
        
#        if (nl==1)
#        {
#            out.all.list <- list(out.all)
#            nstates.list <- list(nstates)
#        }
#        else
#        {
#            out.all.list[[nl]] <- out.all
#            nstates.list[[nl]] <- nstates
#        }
        
        ##cloneinfo <- as.data.frame(cbind(rep(chrom, length(kb.ord)), kb.ord))
        ##dimnames(cloneinfo)[[2]] <- c("Chrom", "kb")
    }
#    list(out.list = out.all.list, nstates.list = nstates.list)
    list(out.list = out.all, nstates.list = nstates)
    
}


#criteria used below are:
	# AIC or BIC
	#if BIC is choosen then it is possible to enter a value for delta as well
"fitHMM" <-
function (MA, vr = 0.01, maxiter = 100, criteria = "AIC", delta = NA, full.output = FALSE) 
{
  #some clunky code so you can put the Criteria argument in characters and still perform the
  #boolean opperators on it below:
          if( criteria == "AIC") {crit = 1}
          else if (criteria == "BIC") {crit = 2}
          else crit = 0
  
    if ((crit == 1) || (crit == 2)) {
    datainfo = MA$genes
    dat = log2.ratios(MA)
    chrom.uniq <- unique(datainfo$Chr)
    nstates <- matrix(NA, nrow = length(chrom.uniq), ncol = ncol(dat))

    #matrix template.  It saves me having to define a new empty matrix 6 times in next bit of code
    	template = matrix(NA,nrow(dat),ncol(dat),dimnames=dimnames(dat))
    #we use template here
    if (full.output == TRUE) {
      segList <- list(M.predicted=template,dispersion=template,state=template,rpred=template,prob=template)
    }
    else {
      segList <- list(M.predicted=template,dispersion=template,state=template,)
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
            cat(j, " ")
            res <- try(states.hmm.func(sample = i, chrom = chrom.uniq[j], 
                dat = dat, datainfo = datainfo, vr = vr, maxiter = maxiter, 
                criteria = crit, delta = delta))
                  nstates[j, i] <- res$nstates.list
 				foo = dat[datainfo$Chr == chrom.uniq[j],i]	
 				segList$M.predicted[(counter+1):(counter+length(foo)),i] = as.matrix(res$out.list[,4])
 				segList$dispersion[(counter+1):(counter+length(foo)),i] = as.matrix(res$out.list[,5])
 		#		segList$obs[(counter+1):(counter+length(foo)),i] = as.matrix(res$out.list[,6])
                                segList$state[(counter+1):(counter+length(foo)),i] = as.matrix(res$out.list[,1])
                                if (full.output == TRUE) { #adding the additional output
                                  
                                  segList$rpred[(counter+1):(counter+length(foo)),i] = as.matrix(res$out.list[,4])
                                  segList$prob[(counter+1):(counter+length(foo)),i] = as.matrix(res$out.list[,5])
                                }
            		counter = counter + length(foo)
        }
        cat("\n")
    }
                segList$M.observed = MA$M
                segList$num.states = nstates
                colnames(segList$num.states) <- colnames(dat)
		segList$genes <- datainfo
		new("SegList",segList)
  }
        else {
          cat("You must enter AIC or BIC for the criteria argument\n")
        }
          
}
