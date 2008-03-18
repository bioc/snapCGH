runTilingArray <- function(input, maxSeg = 5, maxk = 200, criteria = "BIC"){

  #Because both DNAcopy and tilingArray have methods called "segment"
  #I have to detach DNAcopy from the search path in order to run
  #segment for tilingArray.  Ugly, but I can't see another solution at the moment
  if(length(which(search() == "package:DNAcopy")) == 1){
    detach("package:DNAcopy")
  }
  
  chrom.uniq <- unique(input$genes$Chr)

  if(class(input) == "MAList"){
    if (is.null(input$design)) 
      stop("MA$design component is null")

    for(i in 1:length(input$design)){
      temp <- input$design[i]* input$M[,i]
      input$M[,i] <- temp
    }
  }
  
  template <- matrix(NA, nrow(log2ratios(input)), ncol(log2ratios(input)), dimnames = dimnames(log2ratios(input)))

   ### creating the rownames to be used in segList$num.states ####
  rowtemp <- vector()
  rowtemp[1:length(unique(input$genes$Chr))] <- paste("Chrom",unique(input$genes$Chr))
  
  seg.info <- list(M.predicted = template, state = template, M.observed = template,
                   num.states = matrix(NA, length(unique(input$genes$Chr)), ncol(log2ratios(input)), dimnames = list(rowtemp, colnames(input))))

  seg.info$M.observed = input$M
  
  for(i in 1:ncol(input)){
    counter = 0

    for(j in chrom.uniq){

      log2ratios <- input$M[input$genes$Chr == j, i]
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
  seg.info$genes <- input$genes
  seg.info$method <- "Picard"
  #Re-attaching DNAcopy to the search path
  library(DNAcopy, verbose = FALSE, warn.conflicts = FALSE)
  
  new("SegList",seg.info)
}
    
