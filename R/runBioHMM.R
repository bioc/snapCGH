"runBioHMM" <-
function (input, useCloneDists = TRUE, covariates = NULL, criteria="AIC", delta=NA, var.fixed=FALSE, epsilon = 1.0e-6, numiter = 30000) 
{
  if(class(input) == "MAList"){
    if (is.null(input$design)) 
      stop("MA$design component is null")

    for(i in 1:length(input$design)){
      temp <- input$design[i]* input$M[,i]
      input$M[,i] <- temp
    }
  }
  
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

  if ((crit == 1) || (crit == 2)) {
    datainfo = input$genes
    dat = log2ratios(input)   
    chrom.uniq <- unique(datainfo$Chr)
    nstates <- matrix(NA, nrow = length(chrom.uniq), ncol = ncol(dat))

    #matrix template.  It saves me having to define a new empty matrix 6 times in next bit of code
    template = matrix(NA,nrow(dat),ncol(dat),dimnames=dimnames(dat))
    #we use template here
    segList <- list(M.predicted=template,variance=template,state=template)

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
        res <- try(fit.model(sample=i, chrom=chrom.uniq[j], 
                             dat=dat, datainfo=datainfo, useCloneDists = useCloneDists, covariates=covariates, 
               aic = aic, bic = bic, delta=delta, var.fixed=var.fixed, epsilon = epsilon, numiter = numiter))
        nstates[j,i] <- res$nstates.list
        foo = dat[datainfo$Chr == chrom.uniq[j],i]	
        segList$M.predicted[(counter+1):(counter+length(foo)),i] = as.matrix(res$out.list$mean)
        segList$variance[(counter+1):(counter+length(foo)),i] = as.matrix(res$out.list$var)
        segList$state[(counter+1):(counter+length(foo)),i] = as.matrix(res$out.list$state)
        counter = counter + length(foo)
      }
      cat("\n")
    }
    segList$M.observed = input$M
    segList$num.states = nstates
    colnames(segList$num.states) <- colnames(dat)
    rownames(segList$num.states) <- paste("Chrom", unique(input$genes$Chr))
    segList$method <- "BioHMM"
    segList$genes <- datainfo
    new("SegList",segList)
  }
  else {
    cat("You must enter AIC or BIC for the criteria argument\n")
  }
}
