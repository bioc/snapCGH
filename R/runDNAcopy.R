runDNAcopy <- function(MA, smooth.region=2, outlier.SD.scale = 4, smooth.SD.scale = 2, trim=0.025, alpha = 0.01, p.method = c("hybrid", 
    "perm"), kmax = 25, nmin = 200, window.size = NULL, overlap = 0.25, undo.splits = c("none", "prune", "sdundo"), 
    undo.prune = 0.05, undo.SD = 3) {
  
  if (is.null(MA$design)) 
        stop("MA$design component is null")

  for(i in 1:length(MA$design)){
  temp <- MA$design[i]* MA$M[,i]
  MA$M[,i] <- temp
  }

  cna <- CNA(log2ratios(MA), MA$genes$Chr, MA$genes$Position, sampleid = colnames(MA$M.predicted))
  cna <- smooth.CNA(cna, smooth.region=smooth.region, outlier.SD.scale = outlier.SD.scale, smooth.SD.scale = smooth.SD.scale, trim=trim) 
  dna <- segment(cna, alpha = alpha, p.method = p.method, kmax = kmax, nmin = nmin, window.size = window.size, overlap = overlap,
                 trim = trim, undo.splits = undo.splits, undo.prune = undo.prune, undo.SD = undo.SD)
  #changing the output back to the segmentation.info object

  template <- matrix(NA, nrow(log2ratios(MA)), ncol(log2ratios(MA)), dimnames = dimnames(log2ratios(MA)))

   ### creating the rownames to be used in segList$num.states ####
    rowtemp <- vector()
    rowtemp[1:length(unique(MA$genes$Chr))] <- paste("Chrom",unique(MA$genes$Chr))
  
  seg.info <- list(M.predicted = template, state = template, M.observed = template,
                   num.states = matrix(NA, length(unique(MA$genes$Chr)), ncol(log2ratios(MA)), dimnames = list(rowtemp, colnames(MA))))

  nsamples <- ncol(log2ratios(MA))
  names <- unique(dna$output$ID)

  seg.info$M.observed <- log2ratios(MA)
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
  seg.info$genes <- MA$genes
new("SegList",seg.info)
}
