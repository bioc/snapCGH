runGLAD <- function(MA, smoothfunc="aws", base=FALSE, sigma = NULL, bandwidth=10, round=2, lambdabreak=8, lambdacluster=8, lambdaclusterGen=40, type="tricubic", param=c(d=6), alpha=0.001, method="centroid", nmax=8, verbose=FALSE, ...)

{

  if (is.null(MA$design)) 
        stop("MA$design component is null")

  for(i in 1:length(MA$design)){
  temp <- MA$design[i]* MA$object$M[,i]
  MA$M[,i] <- temp
  }
  
  library(GLAD)
    template = matrix(NA,nrow(MA$M),ncol(MA),dimnames=dimnames(MA))
    segList <- list(M.predicted=template,state=template,M.observed=template, num.states=matrix(NA,nrow = length(unique(MA$genes$Chr)), ncol = ncol(MA), dimnames = dimnames(MA)))
    segList$M.observed <- log2.ratios(MA)
    for(k in 1:ncol(log2.ratios(MA))){
      
      cat("Analyzing sample: ", k, "\n")
      profileValues  <- data.frame(PosOrder = (rownames(MA$genes)), LogRatio =(log2.ratios(MA))[,k], PosBase = MA$genes$Position, Chromosome = MA$genes$Chr)
      profileCGH <- list(profileValues=profileValues)
      class(profileCGH) <- "profileCGH"
      
      if (verbose)
      {
        print("GLAD: starting function")
        call <- match.call()
        print(paste("Call function:", call))
      }
    #profileCGH <- list(profileValues=data)	
    #class(profileCGH) <- "profileCGH"	
	
    # Breakpoints detection    
    profileCGH <- chrBreakpoints(profileCGH, smoothfunc=smoothfunc, base=base, sigma=sigma, bandwidth=bandwidth, round=round, verbose=verbose, ...)
    
    # LogRatio are median-centered	
    median <- median(na.omit(profileCGH$profileValues$LogRatio)) 
    profileCGH$profileValues$LogRatio <- profileCGH$profileValues$LogRatio - median	
    profileCGH$profileValues$Smoothing <- profileCGH$profileValues$Smoothing - median	
	    
    profileAux <- NULL	
                                        # profile by chromosome	
    nbzonetot <- 0 #total number of zones that have been previously identify	

    labelChr <- sort(unique(profileCGH$profileValues$Chromosome))	
    for (i in 1:length(labelChr))
      {
        indexChr <- which(profileCGH$profileValues$Chromosome==labelChr[i])	
        subset <- profileCGH$profileValues[indexChr,]

	segList$num.states[i,k] <- length(unique(subset$Smoothing))

        profileChr <- list(profileValues=subset)	
        class(profileChr) <- "profileChr"	
	
        profileChr <- removeBreakpoints(profileChr, lambda=lambdabreak, alpha=alpha, type=type, param=param, verbose=verbose)	
        profileChr <- detectOutliers(profileChr, region="Region", alpha=alpha, verbose=verbose)	
	
	# ça ne doit pas servir : à vérifier
        nmin <- 1	
        if (length(which(profileChr$profileValues$Breakpoints==1))>=1)
          {
            nmin <- 2
          }

        profileChr <- findCluster(profileChr, method=method, genome=FALSE, lambda=lambdacluster, nmin=1, nmax=nmax,type=type, param=param, verbose=verbose)
        profileChr <- detectOutliers(profileChr, region="ZoneChr", alpha=alpha)	
 	
        nbzone <- length(unique((profileChr$profileValues$ZoneChr[which(profileChr$profileValues$ZoneChr!=0)])))	
        profileChr$profileValues$ZoneChr <- profileChr$profileValues$ZoneChr + nbzonetot	
        nbzonetot <- nbzonetot + nbzone	
	
        profileAux <- rbind(profileAux, profileChr$profileValues)
      }
		
    profileCGH$profileValues <- profileAux	
	
    class(profileCGH) <- "profileChr"

    if (verbose) print("GLAD: starting clustering for whole genome")
    profileCGH <- findCluster(profileCGH, region="ZoneChr", method=method, genome=TRUE, lambda=lambdaclusterGen, nmin=1, nmax=nmax, type=type, param=param, verbose=verbose)
    if (verbose)
      {
        print("GLAD: ending clustering for whole genome")
        print("GLAD: starting affectationGNL")      
      }

    class(profileCGH) <- "profileCGH"
    profileCGH <- affectationGNL(profileCGH, verbose=verbose)

    if (verbose)
      {
        print("GLAD: ending affectationGNL")
        print("GLAD: ending function")
      }

#    profileCGH$profileValues <- profileCGH$profileValues[order(profileCGH$profileValues$PosOrder),]
    segList$M.predicted[,k] <- profileCGH$profileValues$Smoothing
    segList$state[,k] <- profileCGH$profileValues$Level
}
    segList$genes <- MA$genes
    new("SegList", segList)
}

