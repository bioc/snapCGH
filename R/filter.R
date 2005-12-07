filterClones <- function(MA, filterFunc, ...){
  remove <- as.matrix(filterFunc(MA, ...))
  for(i in 1:ncol(MA)){
    MA$M[remove[,i],i] <- NA
  }
  MA
}

removeByWeights <- function(MA, weights=MA$weights, threshold = 0.2){

  index <- weights[,] < threshold
  as.matrix(index)
}

