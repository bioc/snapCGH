"convert.output" <-
function(input){
  holder <- list()
  for (i in 1:length(input)){
  holder[[i]] <- list()}
  for(i in 1:length(input)){
    holder[[i]]$genes <- matrix(NA, nrow = length(input[[i]]$clones$mid.point),
                                ncol = 2)
  }
  for(i in 1:length(input)){
    holder[[i]]$M <- as.matrix(input[[i]]$datamatrix)
    holder[[i]]$genes[,1] <- input[[i]]$clones$mid.point
    holder[[i]]$genes[,2] <- rep(input[[i]]$chrom,length(input[[i]]$clones$mid.point))
    
    holder[[i]] <- new("MAList", holder[[i]])
	colnames(holder[[i]]$genes) <- c("Position", "Chr")
    holder[[i]]$genes <- as.data.frame(holder[[1]]$genes)
  }
  holder
}

