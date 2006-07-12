#mod.findBreakPoints <- function(states) 
#{
#    bpoints <- vector()
#    for (i in 2:length(states)) {
#        if (states[i] != states[i - 1]) 
#            bpoints <- c(bpoints, i)
#    }
#    bpoints
#}

#mod1.findBreakPoints <- function(states) 
#{
#    bpoints <- vector()
#    for (i in 2:length(states)) {
#        if (states[i] != states[i - 1]) 
#            bpoints <- c(bpoints,i-1, i, i+1)
#    }
#    bpoints <- bpoints[bpoints <= length(states)]
#    bpoints <- unique(bpoints)
#    bpoints
#}

#mod2.findBreakPoints <- function(states) 
#{
#    bpoints <- vector()
#    for (i in 2:length(states)) {
#        if (states[i] != states[i - 1]) 
#            bpoints <- c(bpoints,i-2,i-1, i, i+1,i+2)
#    }
#    bpoints <- bpoints[bpoints <= length(states)]
#    bpoints <- unique(bpoints)
#    bpoints
#}

compareBreakPoints <- function(states, offset){
  bpoints <- vector()
  for(i in 2:length(states)){
    if(states[i] != states[i-1])
      bpoints <- switch((offset+1), c(bpoints,i), c(bpoints, c(i-1:i)),
                        c(bpoints, c(i-2:i+2)))
  }
  bpoints <- bpoints[bpoints <= length(states)]
  bpoints <- unique(bpoints)
  bpoints
}

# in the function below, TrueSeg is the true Segmented output, 

compareSegmentations <- function(TrueSeg,offset = 0,...){
  
  objects <- list(...)
  nobj <- length(objects)
  nsamps <- ncol(TrueSeg)
  
  
  # The vector TrueOutVec has an entry of 1 if a probe is the first probe of a new state. It takes a value of 0 otherwise.

  TrueOutVec <- vector()
  
  for (i in 1:nsamps){
    pre.temp <- TrueSeg$state[,i]
    for (j in sort(unique(TrueSeg$genes$Chr))){
      temp <- pre.temp[TrueSeg$genes$Chr==j]

      inter <- vector()
      inter <- rep(0,length(temp))
#      BP.func <- switch(offset,
#                        mod.findBreakPoints(temp),
#                        mod1.findBreakPoints(temp),
#                        mod2.findBreakPoints(temp))
      BP.func <- compareBreakPoints(temp, offset)
      inter[BP.func] <- 1    
      TrueOutVec <- c(TrueOutVec,inter)
    }
  }

# segmented values - analogous to TrueOutVec.
  
  SegBP <- list()
  
  for (i in 1:nobj){
    SegBP[[i]] <- vector()

    for (j in 1:nsamps){
      pre.temp <- objects[[i]]$state[,j]
      for (k in sort(unique(objects[[i]]$genes$Chr))){
        temp <- pre.temp[objects[[i]]$genes$Chr==k]
        
        inter <- vector()
        inter <- rep(0,length(temp))
#        BP.func <- switch(offset,
#                        mod.findBreakPoints(temp),
#                        mod1.findBreakPoints(temp),
#                        mod2.findBreakPoints(temp))
        BP.func <- compareBreakPoints(temp, offset)
        inter[BP.func] <- 1    
        SegBP[[i]] <- c(SegBP[[i]],inter)
      }
    }
  }
    
  print("SegOut and TrueOutVec")
  
  
# Finding the TPR and FDR for each genome individually.

  TPR = matrix(ncol = nsamps, nrow = nobj)
  FDR = matrix(ncol = nsamps, nrow = nobj)

  colnames(TPR) <- colnames(FDR) <- colnames(TrueSeg)

  Len <- nrow(TrueSeg$genes)
  
  for (i in 1:nobj){   
    for (j in 1:nsamps){

      inter1 <- SegBP[[i]][(Len*(j-1) + 1) : (Len*j)]
      inter2 <- TrueOutVec[(Len*(j-1) + 1) : (Len*j)]
      inter3.TPR <- ifelse((inter1 ==1) & (inter1 == inter2),1,0)
      inter3.FDR <- ifelse((inter1 ==1) & (inter1 != inter2),1,0) 

      denom.TPR <- sum(inter2)
      denom.FDR <- sum(inter1)
      
      TPR[i,j] <- sum(inter3.TPR)/denom.TPR
      FDR[i,j] <- sum(inter3.FDR)/denom.FDR
  }
  
  print("TPR and FDR calculated")

  }

  output <- list()
  output$TPR <- TPR
  output$FDR <- FDR
  return(output)
}
  