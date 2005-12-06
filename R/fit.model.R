"fit.model" <-
function (sample, chrom, dat, datainfo = clones.info, covariates, vr = 0.01, 
    iterlim = iterlim, criteria = 1, delta = 1, nlists = 1,var.fixed=FALSE) 
{
    obs <- dat[datainfo$Chr == chrom, sample]
    kb <- datainfo$Position[datainfo$Chr == chrom]

    no.cov <- (ncol(covariates)-2)/ncol(dat)
    cov <- covariates[covariates$Chr == chrom,(3 + ((sample-1)*no.cov)):(2 + no.cov*sample)] 
    cov <- as.matrix(cov)
    
    obs.ord <- obs[order(kb)]
    kb.ord <- kb[order(kb)]
    ind.nonna <- which(!is.na(obs.ord))
    y <- obs.ord[ind.nonna] 
    kb <- kb.ord[ind.nonna]
    covars <- cov
    numobs <- length(y)
    
    if(numobs > 5){
     
    k1 <- 2
    k2 <- 8
    k3 <- 15
    k4 <- 24
    k5 <- 35

    z1.pre <- one.state(y,iterlim)
    z2.pre <- two.states(y,covars,iterlim,var.fixed)
    z3.pre <- three.states(y,covars,iterlim,var.fixed)
    z4.pre <- four.states(y,covars,iterlim,var.fixed)
    z5.pre <- five.states(y,covars,iterlim,var.fixed)

    z1 <- find.param.one(z1.pre)
    z2 <- find.param.two(z2.pre,var.fixed)
    z3 <- find.param.three(z3.pre,var.fixed)
    z4 <- find.param.four(z4.pre,var.fixed)
    z5 <- find.param.five(z5.pre,var.fixed)    

    if (length(names(z1)) == 0) {
        z1$minus.logLikelihood <- NA
    }
    if (length(names(z2)) == 0) {
        z2$minus.logLikelihood <- NA
    }
    if (length(names(z3)) == 0) {
        z3$minus.logLikelihood <- NA
    }
    if (length(names(z4)) == 0) {
        z4$minus.logLikelihood <- NA
    }
    if (length(names(z5)) == 0) {
        z5$minus.logLikelihood <- NA
    }

        
    for (nl in 1:nlists) {
        if (criteria == 1) {
            factor <- 2
        }
        else if (criteria == 2) {
          factor <- log(numobs) * delta
        }
        lik <- c((z1$minus.logLikelihood + k1 * factor/2), (z2$minus.logLikelihood + k2 * factor/2), 
            (z3$minus.logLikelihood + k3 * factor/2), (z4$minus.logLikelihood + k4 * factor/2), 
            (z5$minus.logLikelihood + k5 * factor/2))
        likmin <- which.min(lik)
        if (likmin == 1) {
            z <- z1
            name <- "z1"
            nstates <- 1
        }
        else if (likmin == 2) {
            z <- z2
            name <- "z2"
            nstates <- 2
        }
        else if (likmin == 3) {
            z <- z3
            name <- "z3"
            nstates <- 3
        }
        else if (likmin == 4) {
            z <- z4
            name <- "z4"
            nstates <- 4
        }
        else if (likmin == 5) {
            z <- z5
            name <- "z5"
            nstates <- 5
        }
        
        # Fine up to here - this is simply finding the model which best satisfies the optimality criterion chosen
        # The next step is to generate the output from my model. This will require finding the heterogeneous transition matrix.
       
       if (nstates != 1){
         trans.mat <- list()
         for (j in 1:(length(y)-1))
            {trans.mat[[j]] <- z$LH.trans + exp(-(covars[j,1]^(z$rate1))*prod(covars[j,-1]))*z$RH.trans }
          }

        # We then apply the Viterbi algorithm as written by myself to do the rest of the calcualtions. 

        switch(nstates,
               Vit.seg <- rep(1,length(y)),
               Vit.seg <- Viterbi.two(y,z,trans.mat),
               Vit.seg <- Viterbi.three(y,z,trans.mat),
               Vit.seg <- Viterbi.four(y,z,trans.mat),
               Vit.seg <- Viterbi.five(y,z,trans.mat))
        
                
        if (nstates > 1) {
            maxstate.unique <- unique(Vit.seg)
            mean <- rep(0, length(y))
            var <- rep(0, length(y))
            for (m in 1:length(maxstate.unique)) {
                mean[Vit.seg == maxstate.unique[m]] <- mean(y[Vit.seg == maxstate.unique[m]])
                var[Vit.seg == maxstate.unique[m]] <- var(y[Vit.seg == maxstate.unique[m]])
            }
        }
        else {
            maxstate <- rep(1, length(y))
            mean <- rep(mean(y), length(y))
            var <- rep(var(y), length(y))
        }
        out <- cbind(matrix(Vit.seg, ncol = 1), matrix(mean, ncol = 1), matrix(var, ncol = 1))
        out.all <- matrix(NA, nrow = length(kb.ord), ncol = 4)
        out.all[ind.nonna, 1:3] <- out
        out.all[, 4] <- obs.ord
        out.all <- as.data.frame(out.all)
        dimnames(out.all)[[2]] <- c("state", "mean", "var", "obs")
##        if (nl == 1) {
##            out.all.list <- list(out.all)
##            nstates.list <- list(nstates)
##        }
##        else {
##            out.all.list[[nl]] <- out.all
##            nstates.list[[nl]] <- nstates
##        }
      }
  }
   else {
     out.all <- matrix(NA, nrow = length(kb.ord), ncol = 4)
     out.all[ind.nonna,1] <- c(rep(1,numobs))
     out.all[ind.nonna,2] <- c(rep(mean(obs.ord),numobs))
     out.all[ind.nonna,3] <- c(rep(var(obs.ord),numobs))
     out.all[ind.nonna,4] <- obs.ord
     out.all <- as.data.frame(out.all)
     dimnames(out.all)[[2]] <- c("state", "mean", "var", "obs")
     nstates = 1
   }
    list(out.list = out.all, nstates.list = nstates)
  }
