runDNAcopy <- function(input, smooth.region=2, outlier.SD.scale = 4, smooth.SD.scale = 2, trim=0.025, alpha = 0.01, p.method = c("hybrid", 
    "perm"), kmax = 25, nmin = 200, undo.splits = c("none", "prune", "sdundo"), 
    undo.prune = 0.05, undo.SD = 3, nperm=10000, eta=0.05, sbdry=NULL) {

  if(length(which(search() == "package:tilingArray")) == 1){
    detach("package:tilingArray")
  }

  if(class(input) == "MAList"){
    if (is.null(input$design)) 
      stop("MA$design component is null")
    
    for(i in 1:length(input$design)){
      temp <- input$design[i]* input$M[,i]
      input$M[,i] <- temp
    }
  }

  cna <- CNA(log2ratios(input), input$genes$Chr, input$genes$Position, sampleid = colnames(input$M.observed))
  cna <- smooth.CNA(cna, smooth.region=smooth.region, outlier.SD.scale = outlier.SD.scale, smooth.SD.scale = smooth.SD.scale, trim=trim) 
  dna <- segment(cna, alpha = alpha, p.method = p.method, kmax = kmax, nmin = nmin, nperm=nperm, eta = eta, sbdry=sbdry,
                 trim = trim, undo.splits = undo.splits, undo.prune = undo.prune, undo.SD = undo.SD)
  #changing the output back to the segmentation.info object

  template <- matrix(NA, nrow(log2ratios(input)), ncol(log2ratios(input)), dimnames = dimnames(log2ratios(input)))

   ### creating the rownames to be used in segList$num.states ####
    rowtemp <- vector()
    rowtemp[1:length(unique(input$genes$Chr))] <- paste("Chrom",unique(input$genes$Chr))
  
  seg.info <- list(M.predicted = template, state = template, M.observed = template,
                   num.states = matrix(NA, length(unique(input$genes$Chr)), ncol(log2ratios(input)), dimnames = list(rowtemp, colnames(input))))

  nsamples <- ncol(log2ratios(input))
  names <- unique(dna$output$ID)

  seg.info$M.observed <- log2ratios(input)
  for(i in 1:nsamples) {
    temp <- dna$output[dna$output$ID == names[i],]

    counter = 0
    j2 = 0
    for(j in 1:nrow(as.matrix(temp$ID))) {
      nclones <- temp$num.mark[j]
      seg.info$M.predicted[(counter+1):(counter+nclones),i] = temp$seg.mean[j]

      if((j != 1) && (temp$chrom[j] == temp$chrom[j-1])) {
        state = state + 1 
      }
      else {
        state = 1
        j2 = j2+1
      }
      
      seg.info$num.states[j2,i] <- state
      seg.info$state[(counter+1):(counter+nclones),i] = state

      counter = counter + nclones
    }     
  }
  seg.info$method <- "DNAcopy"
  seg.info$genes <- input$genes

    #Re-attaching tilingArray to the search path
  library(tilingArray, verbose = FALSE, warn.conflicts = FALSE)
  
  new("SegList",seg.info)
}
