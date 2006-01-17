"fit.model" <-
function (sample, chrom, dat, datainfo = clones.info, covariates, 
    iterlim = iterlim, aic = TRUE, bic = FALSE, delta = 1, nlists = 1,var.fixed=FALSE) 
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
    data <- y
    
    if(numobs > 5){
    
    k1 <- 2
    k2 <- 8
    k3 <- 15
    k4 <- 24
    k5 <- 35

    one.NM <- function(nrow, xin, data){
	res <- .C("one_state_praxis", as.integer(nrow), as.double(xin), as.double(data), result = double(1), PACKAGE = "snapCGH")
  	list(x = res[[2]],val = res[[4]])
    	}

    two.NM <- function(nrow, xin, data, covars, var.fixed){
      res <- .C("two_states_praxis", as.integer(nrow), as.double(xin), as.double(data), as.double(covars), as.integer(var.fixed), result = double(1), PACKAGE = "snapCGH")
 	list(x = res[[2]],val = res[[6]])
	}

    three.NM <- function(nrow, xin, data, covars, var.fixed){
      res <- .C("three_states_praxis", as.integer(nrow), as.double(xin), as.double(data), as.double(covars), as.integer(var.fixed), result = double(1), PACKAGE = "snapCGH")

	list(x = res[[2]],val = res[[6]])
	}

    four.NM <- function(nrow, xin, data, covars, var.fixed){
      res <- .C("four_states_praxis", as.integer(nrow), as.double(xin), as.double(data), as.double(covars), as.integer(var.fixed), result = double(1), PACKAGE = "snapCGH")
 	list(x = res[[2]],val = res[[6]])
	}

    five.NM <- function(nrow, xin, data, covars, var.fixed){
      res <- .C("five_states_praxis", as.integer(nrow), as.double(xin), as.double(data), as.double(covars), as.integer(var.fixed), result = double(1), PACKAGE = "snapCGH")
 	list(x = res[[2]],val = res[[6]])	
	}

    # Initialisation for one state

    init.mean <- mean(data)

    z1.init <- c(init.mean,0.5)

    # Initialisation for two states

    init.mean.two <- pam(data,2)$medoids

init.var.two <- vector()

if (var.fixed == FALSE){
if (length(pam(data,2)$data[pam(data,2)$clustering==1]) > 1)  init.var.two[1] <- log(sqrt(var(pam(data,2)$data[pam(data,2)$clustering==1]))) else init.var.two[1] <- log(0.5)
if (length(pam(data,2)$data[pam(data,2)$clustering==2]) > 1)  init.var.two[2] <- log(sqrt(var(pam(data,2)$data[pam(data,2)$clustering==2]))) else init.var.two[2] <- log(0.5)} else {
init.var.two[1] <- log(sqrt(var(data)))
init.var.two[2] <- log(sqrt(var(data)))}

    z2.init <- c(init.mean.two[,1],init.var.two,-1,-3.6,-3.6,0)

    # Initialisation for three states

init.mean.three <- vector()
init.mean.three <- pam(data,3)$medoids

init.var.three <- vector()


if (var.fixed == FALSE){
if (length(pam(data,3)$data[pam(data,3)$clustering==1]) > 1)  init.var.three[1] <- log(sqrt(var(pam(data,3)$data[pam(data,3)$clustering==1]))) else init.var.three[1] <- log(0.5)
if (length(pam(data,3)$data[pam(data,3)$clustering==2]) > 1)  init.var.three[2] <- log(sqrt(var(pam(data,3)$data[pam(data,3)$clustering==2]))) else init.var.three[2] <- log(0.5)
if (length(pam(data,3)$data[pam(data,3)$clustering==3]) > 1)  init.var.three[3] <- log(sqrt(var(pam(data,3)$data[pam(data,3)$clustering==3]))) else init.var.three[3] <- log(0.5)} else {
init.var.three[1] <- log(sqrt(var(data)))
init.var.three[2] <- log(sqrt(var(data)))
init.var.three[3] <- log(sqrt(var(data)))}


    z3.init <- c(init.mean.three[,1],init.var.three,-0.7,-0.7,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,0)

    # Initialisation for four states

    init.mean.four <- vector()
init.mean.four <- pam(data,4)$medoids

init.var.four <- list()

init.var.four <- vector()

if (var.fixed == FALSE){
if (length(pam(data,4)$data[pam(data,4)$clustering==1]) > 1)  init.var.four[1] <- log(sqrt(var(pam(data,4)$data[pam(data,4)$clustering==1]))) else init.var.four[1] <- log(0.5)
if (length(pam(data,4)$data[pam(data,4)$clustering==2]) > 1)  init.var.four[2] <- log(sqrt(var(pam(data,4)$data[pam(data,4)$clustering==2]))) else init.var.four[2] <- log(0.5)
if (length(pam(data,4)$data[pam(data,4)$clustering==3]) > 1)  init.var.four[3] <- log(sqrt(var(pam(data,4)$data[pam(data,4)$clustering==3]))) else init.var.four[3] <- log(0.5)
if (length(pam(data,4)$data[pam(data,4)$clustering==4]) > 1)  init.var.four[4] <- log(sqrt(var(pam(data,4)$data[pam(data,4)$clustering==4]))) else init.var.four[4] <- log(0.5)} else {
init.var.four[1] <- log(sqrt(var(data)))
init.var.four[2] <- log(sqrt(var(data)))
init.var.four[3] <- log(sqrt(var(data)))
init.var.four[4] <- log(sqrt(var(data)))}

    z4.init <- c(init.mean.four[,1],init.var.four,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,0)
    
    # Initialisation for five states

    init.mean.five <- vector()
    init.mean.five <- pam(data,5)$medoids

    init.var.five <- list()

    init.var.five <- vector()

    if (var.fixed == FALSE){
      if (length(pam(data,5)$data[pam(data,5)$clustering==1]) > 1)  init.var.five[1] <- log(sqrt(var(pam(data,5)$data[pam(data,5)$clustering==1]))) else init.var.five[1] <- log(0.5)
    if (length(pam(data,5)$data[pam(data,5)$clustering==2]) > 1)  init.var.five[2] <- log(sqrt(var(pam(data,5)$data[pam(data,5)$clustering==2]))) else init.var.five[2] <- log(0.5)
    if (length(pam(data,5)$data[pam(data,5)$clustering==3]) > 1)  init.var.five[3] <- log(sqrt(var(pam(data,5)$data[pam(data,5)$clustering==3]))) else init.var.five[3] <- log(0.5)
    if (length(pam(data,5)$data[pam(data,5)$clustering==4]) > 1)  init.var.five[4] <- log(sqrt(var(pam(data,5)$data[pam(data,5)$clustering==4]))) else init.var.five[4] <- log(0.5)
    if (length(pam(data,5)$data[pam(data,5)$clustering==5]) > 1)  init.var.five[5] <- log(sqrt(var(pam(data,5)$data[pam(data,5)$clustering==5]))) else init.var.five[5] <- log(0.5)} else {
      init.var.five[1] <- log(sqrt(var(data)))
      init.var.five[2] <- log(sqrt(var(data)))
      init.var.five[3] <- log(sqrt(var(data)))
      init.var.five[4] <- log(sqrt(var(data)))
      init.var.five[5] <- log(sqrt(var(data)))}

z5.init <- c(init.mean.five[,1],init.var.five,-0.7,-0.7,-0.7,-0.7,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,0) 
    

    z1.pre <- one.NM(numobs,z1.init ,y)
    z2.pre <- two.NM(numobs, z2.init, y, covars, var.fixed)
    z3.pre <- three.NM(numobs,z3.init,y,covars, var.fixed)
    z4.pre <- four.NM(numobs,z4.init,y,covars, var.fixed)
    z5.pre <- five.NM(numobs,z5.init,y,covars, var.fixed)
    
    #z1.pre <- one.state(y,iterlim)
    #z2.pre <- two.states(y,covars,iterlim,var.fixed)
    #z3.pre <- three.states(y,covars,iterlim,var.fixed)
    #z4.pre <- four.states(y,covars,iterlim,var.fixed)
    #z5.pre <- five.states(y,covars,iterlim,var.fixed)

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
        if ((aic) && (nl == 1)) {
            factor <- 2
        }
        else if (bic) {
            if (aic) {
                factor <- log(numobs) * delta[nl - 1]
            }
            else {
                factor <- log(numobs) * delta[nl]
            }
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
#        if (nl == 1) {
#            out.all.list <- list(out.all)
#            nstates.list <- list(nstates)
#        }
#        else {
#            out.all.list[[nl]] <- out.all
#            nstates.list[[nl]] <- nstates
#        }
      }
  }  else {
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


