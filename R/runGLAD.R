runGLAD <- function(MA, smoothfunc="aws", base=FALSE, sigma = NULL, bandwidth=10, round=2, lambdabreak=8, lambdacluster=8, lambdaclusterGen=40, type="tricubic", param=c(d=6), alpha=0.001, method="centroid", nmax=8, verbose=FALSE, ...)

{

  if (is.null(MA$design)) 
        stop("MA$design component is null")

  for(i in 1:length(MA$design)){
  temp <- MA$design[i]* MA$M[,i]
  MA$M[,i] <- temp
  }
  
    template = matrix(NA,nrow(MA$M),ncol(MA),dimnames=dimnames(MA))
    segList <- list(M.predicted=template,state=template,M.observed=template, num.states=matrix(NA,nrow = length(unique(MA$genes$Chr)), ncol = ncol(MA), dimnames = dimnames(MA)))
    segList$M.observed <- log2ratios(MA)
    for(k in 1:ncol(log2ratios(MA))){
      
      cat("Analyzing sample: ", k, "\n")
      profileValues  <- data.frame(PosOrder = 1:nrow(MA$genes), LogRatio =(log2ratios(MA))[,k], PosBase = MA$genes$Position, Chromosome = MA$genes$Chr)
      profileCGH <- list(profileValues=profileValues)
      class(profileCGH) <- "profileCGH"
      
    out <- glad(profileCGH, smoothfunc=smoothfunc, base=base, sigma=sigma, bandwidth=bandwidth, round=round, lambdabreak=lambdabreak, lambdacluster=lambdacluster, lambdaclusterGen=lambdaclusterGen, type=type, param=param, alpha=alpha, method=method, nmax=nmax, verbose=verbose)
    
#    profileCGH$profileValues <- profileCGH$profileValues[order(profileCGH$profileValues$PosOrder),]
    segList$M.predicted[,k] <- out$profileValues$Smoothing
    segList$state[,k] <- out$profileValues$Level
      for(i in 1:length(unique(MA$genes$Chr))){
        segList$num.states[i,k] <- length(unique(out$profileValues$Level[MA$genes$Chr == i]))
      }
}
    segList$genes <- MA$genes
    new("SegList", segList)

}
