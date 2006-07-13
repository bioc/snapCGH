runPicard <- function(MA, maxSeg, maxk = 200, criteria = "BIC"){

  #Because both DNAcopy and tilingArray have methods called "segment"
  #I have to detach DNAcopy from the search path in order to run
  #segment for tilingArray.  Ugly, but I can't see another solution at the moment
  if(length(which(search() == "package:DNAcopy")) == 1){
    detach("package:DNAcopy")
  }
  
  chrom.uniq <- unique(MA$genes$Chr)

  if (is.null(MA$design)) 
    stop("MA$design component is null")

  for(i in 1:length(MA$design)){
    temp <- MA$design[i]* MA$M[,i]
    MA$M[,i] <- temp
  }
  
  template <- matrix(NA, nrow(log2ratios(MA)), ncol(log2ratios(MA)), dimnames = dimnames(log2ratios(MA)))

   ### creating the rownames to be used in segList$num.states ####
  rowtemp <- vector()
  rowtemp[1:length(unique(MA$genes$Chr))] <- paste("Chrom",unique(MA$genes$Chr))
  
  seg.info <- list(M.predicted = template, state = template, M.observed = template,
                   num.states = matrix(NA, length(unique(MA$genes$Chr)), ncol(log2ratios(MA)), dimnames = list(rowtemp, colnames(MA))))

  seg.info$M.observed = MA$M
  
  for(i in 1:ncol(MA)){
    counter = 0

    for(j in chrom.uniq){

      log2ratios <- MA$M[MA$genes$Chr == j, i]
      seg <- segment(log2ratios, maxk = length(log2ratios), maxseg = min(length(log2ratios), maxSeg))
      selected <- which.max(logLik(seg, penalty = criteria))

      est <- c(1,seg@breakpoints[[selected]][,1], (length(log2ratios)+1))
      counter2 = 0
      statecounter = 1
      for(k in 2:length(est)){
        seg.info$M.predicted[(counter+est[k-1]):(counter+est[k]-1),i] <- mean(log2ratios[(est[k-1]):est[k]-1])
        seg.info$state[(counter+est[k-1]):(counter+est[k]-1),i] <- statecounter
        counter2 = est[k]-1
        statecounter = statecounter + 1
      }
      counter = counter + (est[length(est)]-1)
      seg.info$num.states[j,i] <- (length(est)-1)
    }
  }
  seg.info$genes <- MA$genes
  seg.info$method <- "Picard"
  #Re-attaching DNAcopy to the search path
  library(DNAcopy, verbose = FALSE, warn.conflicts = FALSE)
  
  new("SegList",seg.info)
}
    
