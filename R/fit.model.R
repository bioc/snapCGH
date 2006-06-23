"fit.model" <-
function (sample, chrom, dat, datainfo = clones.info, covariates = NULL, 
          aic = TRUE, bic = FALSE, delta = 1, 
          var.fixed=FALSE, epsilon = 1.0e-6, numiter = 30000) {

    obs <- dat[datainfo$Chr == chrom, sample]
    kb <- datainfo$Position[datainfo$Chr == chrom]

   ##extracting distances from the $genes matrix
    dists.pre = kb[2:length(kb)] - kb[1:(length(kb)-1)]
    dists = dists.pre/(max(dists.pre))
    covars <- as.matrix(dists)
   
    if(!is.null(covariates)){
      no.cov <- ((ncol(covariates)-2)/ncol(dat))
      cov <- covariates[covariates$Chr == chrom,(3 + ((sample-1)*no.cov)):(2 + no.cov*sample)]
      covars <- as.matrix(cbind(covars, cov))
    }
     
#   covars <- as.matrix(cbind(dists, cov))
    
    obs.ord <- obs[order(kb)]
    kb.ord <- kb[order(kb)]
    ind.nonna <- which(!is.na(obs.ord))
    data <- obs.ord[ind.nonna] 
    kb <- kb.ord[ind.nonna]
#    covars <- cov
    numobs <- length(data)
#    data <- y
    
    if(numobs > 5){

      if(var.fixed == FALSE){
        k1 <- 2
        k2 <- 8
        k3 <- 15
        k4 <- 24
        k5 <- 35
      }
      else{
        k1 <- 2
        k2 <- 7
        k3 <- 13
        k4 <- 21
        k5 <- 31
      }
    
#    one.NM <- function(nrow, xin, data, maxfunc, tol){
#	res <- .C("one_state_praxis", as.integer(nrow), as.double(xin), as.double(data), result = double(1), as.integer(maxfunc), as.double(tol), PACKAGE = "snapCGH")
#  	list(x = res[[2]],val = res[[4]])
 #   	}

#    two.NM <- function(nrow, xin, data, covars, var.fixed, epsilon){
#      res <- .C("two_states_nelder", as.integer(nrow), as.double(xin), xout = double(8), result = double(1), as.double(data), as.double(covars), as.integer(var.fixed),  as.double(epsilon), as.integer(0), as.integer(numiter), PACKAGE = "altsnapCGH")
# 	list(x = res[[3]],val = res[[4]])
#	}

#    three.NM <- function(nrow, xin, data, covars, var.fixed, maxfunc, tol){
#      res <- .C("three_states_nelder", as.integer(nrow), as.double(xin), xout = double(15), result = double(1), as.double(data), as.double(covars), as.integer(var.fixed),  as.double(epsilon), as.integer(0), as.integer(numiter), PACKAGE = "altsnapCGH")
#	list(x = res[[3]],val = res[[4]])
#	}

#    four.NM <- function(nrow, xin, data, covars, var.fixed, maxfunc, tol){
#      res <- .C("four_states_nelder", as.integer(nrow), as.double(xin), xout = double(24), result = double(1), as.double(data), as.double(covars), as.integer(var.fixed),  as.double(epsilon), as.integer(0), as.integer(numiter), PACKAGE = "altsnapCGH")
# 	list(x = res[[3]],val = res[[4]])
#	}

#    five.NM <- function(nrow, xin, data, covars, var.fixed, maxfunc, tol){
#      res <- .C("five_states_nelder", as.integer(nrow), as.double(xin), xout = double(35), result = double(1), as.double(data), as.double(covars), as.integer(var.fixed),  as.double(epsilon), as.integer(0), as.integer(numiter), PACKAGE = "altsnapCGH")
# 	list(x = res[[3]],val = res[[4]])	
#	}


    ######### Initialisation for one state ########
    
    init.mean <- mean(data)

    z1.init <- c(init.mean,log(sqrt(var(data))))

    ######## Initialisation for two states ########

      temp2 <- clara(data,2)

    init.mean.two <- temp2$medoids

    init.var.two <- vector()

    if (var.fixed == FALSE){
      for(i in 1:2){
        if (length(temp2$data[temp2$clustering==i]) > 1)
          init.var.two[i] <- log(sqrt(var(temp2$data[temp2$clustering==i])))
        else init.var.two[i] <- log(0.5)
      }
    }
    else{
      init.var.two[1:2] <- log(sqrt(var(data)))
    }

    z2.init <- c(init.mean.two[,1],init.var.two,-1,-3.6,-3.6,0)

    ########## Initialisation for three states ########

      temp3 <- clara(data,3)
      
    init.mean.three <- vector(length = 3)
    init.mean.three <- temp3$medoids

    init.var.three <- vector(length = 3)

    if (var.fixed == FALSE){
      for(i in 1:3){
        if (length(temp3$data[temp3$clustering==i]) > 1)
          init.var.three[i] <- log(sqrt(var(temp3$data[temp3$clustering==i])))
        else init.var.three[i] <- log(0.5)
      }
    }
    else{
      init.var.three[1:3] <- log(sqrt(var(data)))
    }

    z3.init <- c(init.mean.three[,1],init.var.three,-0.7,-0.7,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,0)

    ######## Initialisation for four states #########

      temp4 <- clara(data,4)
      
    init.mean.four <- vector()
    init.mean.four <- temp4$medoids

    init.var.four <- vector()

    if (var.fixed == FALSE){
      for(i in 1:4){
        if (length(temp4$data[temp4$clustering==i]) > 1)
          init.var.four[i] <- log(sqrt(var(temp4$data[temp4$clustering==i])))
        else init.var.four[i] <- log(0.5)
      }
    }
    else {
      init.var.four[1:4] <- log(sqrt(var(data)))
    }
    
    z4.init <- c(init.mean.four[,1],init.var.four,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,0)
    
    ########## Initialisation for five states #############

      temp5 <- clara(data, 5)
      
    init.mean.five <- vector()
    init.mean.five <- temp5$medoids

    init.var.five <- vector()

    if (var.fixed == FALSE){
      for(i in 1:5){
        if(length(temp5$data[temp5$clustering==i]) > 1)
          init.var.five[i] <- log(sqrt(var(temp5$data[temp5$clustering==i])))
        else init.var.five[i] <- log(0.5)
      }
    }
    else {
      init.var.five[1:5] <- log(sqrt(var(data)))
    }

    z5.init <- c(init.mean.five[,1],init.var.five,-0.7,-0.7,-0.7,-0.7,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,0) 

#    z1.pre <- list(x = vector())
#      z1.pre$x[1] <- mean(data)
#      z1.pre$x[2] <- (1/length(data))*(sum((data - mean(data))^2))
#      z1.pre$val <- -sum(log(dnorm(data,z1.pre$x[1],sqrt(z1.pre$x[2]),0)))

#    z2.pre <- two.NM(numobs, z2.init, data, covars, var.fixed, epsilon)
#      print(z2.pre)
#    z3.pre <- three.NM(numobs,z3.init,data,covars, var.fixed, epsilon)
                                        #      print(z3.pre)
#    z4.pre <- four.NM(numobs,z4.init,data,covars, var.fixed, epsilon)
#      print(z4.pre)
#    z5.pre <- five.NM(numobs,z5.init,data,covars, var.fixed, epsilon)
#      print(z5.pre)

      z.pre <- vector("list", length = 5)

      z.pre[[1]]$x <- vector()
      z.pre[[1]]$x[1] <- mean(data)
      z.pre[[1]]$x[2] <- (1/length(data))*(sum((data - mean(data))^2))

    z.pre[[1]]$val <- -sum(log(dnorm(data,z.pre[[1]]$x[1],sqrt(z.pre[[1]]$x[2]),0)))

#    print(z.pre[[1]])
    
    for(i in 2:5){
      init <- paste("z",i,".init", sep="")
      z.pre[[i]] <- run.nelder(numobs, get(init), data, covars, var.fixed, epsilon, numiter, i)
#	print(z.pre[[i]])
    } 
    
    z1 <- find.param.one(z.pre[[1]])
    if (!is.nan(z.pre[[2]]$x[1]))
      z2 <- find.param.two(z.pre[[2]],var.fixed)
    else {
      z2 <- NULL
    }
    if (!is.nan(z.pre[[3]]$x[1]))
      z3 <- find.param.three(z.pre[[3]],var.fixed)
    else {
      z3 <- NULL
    }
    if (!is.nan(z.pre[[4]]$x[1]))
      z4 <- find.param.four(z.pre[[4]],var.fixed)
    else {
      z4 <- NULL
    }
    if (!is.nan(z.pre[[5]]$x[1]))
      z5 <- find.param.five(z.pre[[5]],var.fixed)
    else {
      z5 <- NULL
    }

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
        
    if (aic) {
      factor <- 2
    }
    else if (bic) {
      factor <- log(numobs) * delta
    }
    else {
      stop("No criteria selected")
    }
    
        lik <- c((z1$minus.logLikelihood + k1 * factor/2), (z2$minus.logLikelihood + k2 * factor/2), 
            (z3$minus.logLikelihood + k3 * factor/2), (z4$minus.logLikelihood + k4 * factor/2), 
            (z5$minus.logLikelihood + k5 * factor/2))
      
#      cat("Lik: ",lik, "\n")
        likmin <- which.min(lik)
#      cat("Likmin:", likmin, "\n")
      
        switch(likmin, z <- z1,z <- z2, z <- z3, z <- z4, z <- z5)
        nstates <- likmin
      
        # Fine up to here - this is simply finding the model which best satisfies the optimality criterion chosen
        # The next step is to generate the output from my model. This will require finding the heterogeneous transition matrix.
       
       if (nstates != 1){
         trans.mat <- list()
         for (j in 1:(length(data)-1))
            {trans.mat[[j]] <- z$LH.trans + exp(-(covars[j,1]^(z$rate1))*prod(covars[j,-1]))*z$RH.trans }
          }

        # We then apply the Viterbi algorithm as written by myself to do the rest of the calcualtions. 

        switch(nstates,
               Vit.seg <- rep(1,length(data)),
               Vit.seg <- Viterbi.two(data,z,trans.mat),
               Vit.seg <- Viterbi.three(data,z,trans.mat),
               Vit.seg <- Viterbi.four(data,z,trans.mat),
               Vit.seg <- Viterbi.five(data,z,trans.mat))        
                
        if (nstates > 1) {
            maxstate.unique <- unique(Vit.seg)
            mean <- rep(0, length(data))
            var <- rep(0, length(data))
            for (m in 1:length(maxstate.unique)) {
                mean[Vit.seg == maxstate.unique[m]] <- mean(data[Vit.seg == maxstate.unique[m]])
                var[Vit.seg == maxstate.unique[m]] <- var(data[Vit.seg == maxstate.unique[m]])
            }
        }
        else {
            maxstate <- rep(1, length(data))
            mean <- rep(mean(data), length(data))
            var <- rep(var(data), length(data))
        }
        out <- cbind(matrix(Vit.seg, ncol = 1), matrix(mean, ncol = 1), matrix(var, ncol = 1))
        out.all <- matrix(NA, nrow = length(kb.ord), ncol = 4)
        out.all[ind.nonna, 1:3] <- out
        out.all[, 4] <- obs.ord
        out.all <- as.data.frame(out.all)
        dimnames(out.all)[[2]] <- c("state", "mean", "var", "obs")
        numstates <- length(unique(Vit.seg))
    }
       else {
         out.all <- matrix(NA, nrow = length(kb.ord), ncol = 4)
         out.all[ind.nonna,1] <- c(rep(1,numobs))
         out.all[ind.nonna,2] <- c(rep(mean(obs.ord),numobs))
         out.all[ind.nonna,3] <- c(rep(var(obs.ord),numobs))
         out.all[ind.nonna,4] <- obs.ord
         out.all <- as.data.frame(out.all)
         dimnames(out.all)[[2]] <- c("state", "mean", "var", "obs")
         numstates = 1
       }
    list(out.list = out.all, nstates.list = numstates)
  }



run.nelder <- function(nrow, xin, data, covars, var.fixed, epsilon, numiter,  nstates){
  xout = double(length(xin))
  if(ncol(covars) == 1){
    covars1 = covars[,1]
    covars2 = vector(length = nrow(covars))
    covars3 = vector(length = nrow(covars))
    ncovars = 1
  }
  else if(ncol(covars) == 2){
    covars1 = covars[,1]
    covars2 = covars[,2]
    covars3 = vector(length = nrow(covars))
    ncovars = 2
  }
  else if(ncol(covars) > 2){
    covars1 = covars[,1]
    covars2 = covars[,2]
    covars3 = covars[,3]
    ncovars = 3
  }
  
  res <- .C("runNelderMead", as.integer(nrow), as.double(xin), as.double(xout), result = double(1), as.double(data), as.double(covars1), as.double(covars2), as.double(covars3), as.integer(ncovars), as.integer(var.fixed),  as.double(epsilon), trace = as.integer(0), as.integer(numiter), as.integer(nstates), PACKAGE = "snapCGH")
  list(x = res[[3]],val = res[[4]])
}
