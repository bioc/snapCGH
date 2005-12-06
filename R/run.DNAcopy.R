run.DNAcopy <- function(MA) {
##  library(DNAcopy)
  cna <- CNA(log2.ratios(MA), MA$genes$Chr, MA$genes$Position, sampleid = colnames(MA$M))
  cna <- smooth.CNA(cna) #add some param's later
  dna <- segment(cna)
  #changing the output back to the segmentation.info object

  template <- matrix(NA, nrow(log2.ratios(MA)), ncol(log2.ratios(MA)), dimnames = dimnames(log2.ratios(MA)))
  
  seg.info <- list(M.predicted = template, state = template, M.observed = template, num.states = matrix(NA, length(unique(MA$genes$Chr)), ncol(log2.ratios(MA))))

  nsamples <- ncol(log2.ratios(MA))
  names <- unique(dna$output$ID)

  seg.info$M.observed <- log2.ratios(MA)
  for(i in 1:nsamples) {
    temp <- dna$output[dna$output$ID == names[i],]

    counter = 0
    for(j in 1:nrow(as.matrix(temp$ID))) {
      nclones <- temp$num.mark[j]
      seg.info$M.predicted[(counter+1):(counter+nclones),i] = temp$seg.mean[j]

      if((j != 1) && (temp$chrom[j] == temp$chrom[j-1])) {
        state = state + 1 
      }
      else {
        state = 1
      }
      
      seg.info$num.states[temp$chrom[j],i] <- state
      seg.info$state[(counter+1):(counter+nclones),i] = state

      counter = counter + nclones
    }     
  }
  seg.info$genes <- MA$genes
new("SegList",seg.info)
}
