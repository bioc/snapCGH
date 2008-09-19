runGLAD <- function(input, smoothfunc="lawsglad", base=FALSE, sigma = NULL, bandwidth=10, round=2, lambdabreak=8, lambdacluster=8, lambdaclusterGen=40, type="tricubic", param=c(d=6), alpha=0.001, method="centroid", nmax=8, verbose=FALSE, ...){

  if(class(input) == "MAList"){
    if (is.null(input$design)) 
      stop("MA$design component is null")

    for(i in 1:length(input$design)){
      temp <- input$design[i]* input$M.observed[,i]
      input$M.observed[,i] <- temp
    }
  }
  
 ### creating the rownames to be used in segList$num.states ####
    rowtemp <- vector()
    rowtemp[1:length(unique(input$genes$Chr))] <- paste("Chrom",unique(input$genes$Chr))
  
    template = matrix(NA,nrow(input$M),ncol(input$M),dimnames=dimnames(input$M))
    segList <- list(M.predicted=template,state=template,M.observed=template,
                    num.states=matrix(NA, length(unique(input$genes$Chr)),
                      ncol = ncol(log2ratios(input)), dimnames = list(rowtemp, colnames(input))))
    segList$M.observed <- log2ratios(input)

    for(k in 1:ncol(log2ratios(input))){
      
      cat("Analyzing sample: ", k, "\n")
      profileValues  <- data.frame(PosOrder = 1:nrow(input$genes), LogRatio =(log2ratios(input))[,k], PosBase = input$genes$Position, Chromosome = input$genes$Chr)
      profileCGH <- list(profileValues=profileValues)
      class(profileCGH) <- "profileCGH"
      
    out <- glad(profileCGH, smoothfunc=smoothfunc, base=base, sigma=sigma, bandwidth=bandwidth, round=round, lambdabreak=lambdabreak, lambdacluster=lambdacluster, lambdaclusterGen=lambdaclusterGen, type=type, param=param, alpha=alpha, method=method, nmax=nmax, verbose=verbose)

    segList$M.predicted[,k] <- out$profileValues$Smoothing
    segList$state[,k] <- out$profileValues$Level
      for(i in 1:length(unique(input$genes$Chr))){

        segList$num.states[i,k] <- length(unique(out$profileValues$Level[input$genes$Chr == i]))
      }
    }
  segList$method <- "GLAD"
  segList$genes <- input$genes
  new("SegList", segList)
}
