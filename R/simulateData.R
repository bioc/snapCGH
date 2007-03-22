"simulateData" <-
function(nArrays = 20,
         chrominfo = NULL,
         prb.short.tiled = 0.5,
         prb.long.tiled = 0.5,
         non.tiled.lower.res = 0.9,
         non.tiled.upper.res = 1.1,
         length.clone.lower = 0.05,
         length.clone.upper = 0.2,
         tiled.lower.res = -0.05,
         tiled.upper.res = 0,
         sd = NULL,
         output = FALSE,
         prb.proportion.tiled = c(0.2, 0.2, 0.2, 0.2, 0.2),
         zerolengthnontiled = NULL,
         zerolengthtiled = NULL,
         nonzerolengthnontiled = NULL,
         nonzerolengthtiled = NULL,
         seed = 1){
  
# Randomly select one of the 22 autosomal chromsomes to be our simulation template.
  
  set.seed(seed)

  short.tiled <- vector()
  long.tiled <- vector()
  length.short.tiled <- vector()
  length.long.tiled <- vector()
  start.short.tiled <- vector()
  start.long.tiled <- vector()
  length.chrom.short <- vector()
  length.chrom.long <- vector()

  for (i in 1:22){
    short.tiled[i] <- sample(c(0,1),1,prob=c(1-prb.short.tiled,prb.short.tiled))
    long.tiled[i] <- sample(c(0,1),1,prob=c(1-prb.long.tiled,prb.long.tiled))
    length.chrom.short[i] <- chrominfo.Mb$centromere[i]
    length.chrom.long[i] <- (chrominfo.Mb$length[i] - chrominfo.Mb$centromere[i])
 
    if (short.tiled[i] == 1){
      length.short.tiled[i] <- sample(c(0.2,0.3,0.4,0.5,1),1,prob=prb.proportion.tiled)*length.chrom.short[i]
      if ((length.short.tiled[i]/length.chrom.short[i]) == 1)
        start.short.tiled[i] <- 0
      else 
        start.short.tiled[i] <- 2 + runif(1,0.2,1)*(length.chrom.short[i]-length.short.tiled[i])
    }
    else {
      length.short.tiled[i] <- NA
      start.short.tiled[i] <- NA
    }

    if (long.tiled[i] == 1){
      length.long.tiled[i] <- sample(c(0.2,0.3,0.4,0.5,1),1,prob=prb.proportion.tiled)*length.chrom.long[i]
      if ((length.long.tiled[i]/length.chrom.long[i]) == 1){
        start.long.tiled[i] <- 0}
      else {
        start.long.tiled[i] <- 2 + runif(1,0.2,1)*(length.chrom.long[i]-length.long.tiled[i])
      }
    }
    else {
      length.long.tiled[i] <- NA
      start.long.tiled[i] <- NA
    }
  }
     
  start.pos.short <- list()
  end.pos.short <- list()

  for (i in 1:22){
    start.pos.short[[i]] <- vector()
    end.pos.short[[i]] <- vector()
    start.pos.short[[i]][1] <- 1
    end.pos.short[[i]][1] <- 1.1
  }

  for (i in 1:22){
    if (short.tiled[i] == 0){
      j <- 1
      repeat {
      if (end.pos.short[[i]][j] >= length.chrom.short[i]) {break}    
      j <- j + 1
      start.pos.short[[i]][j] <- end.pos.short[[i]][j-1] + round(runif(1,non.tiled.lower.res,non.tiled.upper.res),4)
      end.pos.short[[i]][j] <- start.pos.short[[i]][j] + round(runif(1,length.clone.lower,length.clone.upper),4)
      }
      start.pos.short[[i]] <- start.pos.short[[i]][1:j-1]
      end.pos.short[[i]] <- end.pos.short[[i]][1:j-1]
    }
  }

  for (i in 1:22){
    if (short.tiled[i] == 1){
      if (length.short.tiled[i] != length.chrom.short[i]){
      j <- 1
      repeat {
        if (end.pos.short[[i]][j] >= start.short.tiled[i]) {break}
        j <- j + 1
        start.pos.short[[i]][j] <- end.pos.short[[i]][j-1] + round(runif(1,non.tiled.lower.res,non.tiled.upper.res),4)
        end.pos.short[[i]][j] <- start.pos.short[[i]][j] + round(runif(1,length.clone.lower,length.clone.upper),4)
        }
        start.pos.short[[i]] <- start.pos.short[[i]][1:j-1]
        end.pos.short[[i]] <- end.pos.short[[i]][1:j-1]
      j <- length(end.pos.short[[i]])
      repeat {
         if (end.pos.short[[i]][j] >= (start.short.tiled[i] + length.short.tiled[i])) {break}
         j <- j + 1
         start.pos.short[[i]][j] <- end.pos.short[[i]][j-1] + round(runif(1,tiled.lower.res,tiled.upper.res),4)
         end.pos.short[[i]][j] <- start.pos.short[[i]][j] + round(runif(1,length.clone.lower,length.clone.upper),4)
         }
         start.pos.short[[i]] <- start.pos.short[[i]][1:j-1]
         end.pos.short[[i]] <- end.pos.short[[i]][1:j-1]
      j <- length(end.pos.short[[i]])
      repeat {
          if (end.pos.short[[i]][j] >= length.chrom.short[i]) break
          j <- j + 1
          start.pos.short[[i]][j] <- end.pos.short[[i]][j-1] + round(runif(1,non.tiled.lower.res,non.tiled.upper.res),4)
          end.pos.short[[i]][j] <- start.pos.short[[i]][j] + round(runif(1,length.clone.lower,length.clone.upper),4)
          }
          start.pos.short[[i]] <- start.pos.short[[i]][1:j-1]
          end.pos.short[[i]] <- end.pos.short[[i]][1:j-1]
    }}}

for (i in 1:22){
  if (short.tiled[i] == 1){
    if (length.short.tiled[i] == length.chrom.short[i]){
      j <- 1
      repeat{
          if (end.pos.short[[i]][j] >= length.chrom.short[i]) {break}
          j <- j + 1
          start.pos.short[[i]][j] <- end.pos.short[[i]][j-1] + round(runif(1,tiled.lower.res,tiled.upper.res),4)
          end.pos.short[[i]][j] <- start.pos.short[[i]][j] + round(runif(1,length.clone.lower,length.clone.upper),4)
          }
       start.pos.short[[i]] <- start.pos.short[[i]][1:j-1]
       end.pos.short[[i]] <- end.pos.short[[i]][1:j-1]
  }}}

# Long arm

  start.pos.long <- list()
  end.pos.long <- list()
  
  for (i in 1:22){
    start.pos.long[[i]] <- vector()
    end.pos.long[[i]] <- vector()
    start.pos.long[[i]][1] <- chrominfo.Mb$centromere[i] + 1
    end.pos.long[[i]][1] <- chrominfo.Mb$centromere[i] + 1.1
  }
  
  for (i in 1:22){
    if (long.tiled[i] == 0){
      j <- 1
      repeat {
      if (end.pos.long[[i]][j] >= (length.chrom.short[i] + length.chrom.long[i])) {break}  
      j <- j + 1
      start.pos.long[[i]][j] <- end.pos.long[[i]][j-1] + round(runif(1,non.tiled.lower.res,non.tiled.upper.res),4)
      end.pos.long[[i]][j] <- start.pos.long[[i]][j] + round(runif(1,length.clone.lower,length.clone.upper),4)
      }
      start.pos.long[[i]] <- start.pos.long[[i]][1:j-1]
      end.pos.long[[i]] <- end.pos.long[[i]][1:j-1]
    }}

for (i in 1:22){
    if (long.tiled[i] == 1){
      if (length.long.tiled[i] != length.chrom.long[i]){
      j <- 1
      repeat {
        if (end.pos.long[[i]][j] >= (length.chrom.short[i] + start.long.tiled[i])) {break}
        j <- j + 1
        start.pos.long[[i]][j] <- end.pos.long[[i]][j-1] + round(runif(1,non.tiled.lower.res,non.tiled.upper.res),4)
        end.pos.long[[i]][j] <- start.pos.long[[i]][j] + round(runif(1,length.clone.lower,length.clone.upper),4)
        }
        start.pos.long[[i]] <- start.pos.long[[i]][1:j-1] 
        end.pos.long[[i]] <- end.pos.long[[i]][1:j-1]
        j <- length(end.pos.long[[i]])
      repeat {
         if (end.pos.long[[i]][j] >= (length.chrom.short[i] + start.long.tiled[i] + length.long.tiled[i])) {break}
         j <- j + 1
         start.pos.long[[i]][j] <- end.pos.long[[i]][j-1] + round(runif(1,tiled.lower.res,tiled.upper.res),4)
         end.pos.long[[i]][j] <- start.pos.long[[i]][j] + round(runif(1,length.clone.lower,length.clone.upper),4)
         }
         start.pos.long[[i]] <- start.pos.long[[i]][1:j-1]
         end.pos.long[[i]] <- end.pos.long[[i]][1:j-1]
         j <- length(end.pos.long[[i]])
      repeat {
          if (end.pos.long[[i]][j] >= (length.chrom.short[i] + length.chrom.long[i])) {break}
          j <- j + 1
          start.pos.long[[i]][j] <- end.pos.long[[i]][j-1] + round(runif(1,non.tiled.lower.res,non.tiled.upper.res),4)
          end.pos.long[[i]][j] <- start.pos.long[[i]][j] + round(runif(1,length.clone.lower,length.clone.upper),4)
          }
          start.pos.long[[i]] <- start.pos.long[[i]][1:j-1] 
          end.pos.long[[i]] <- end.pos.long[[i]][1:j-1]
    }}}

for (i in 1:22){
  if (long.tiled[i] == 1){
    if (length.long.tiled[i] == length.chrom.long[i]){
      j <- 1
      repeat{
          if (end.pos.long[[i]][j] >= (length.chrom.short[i] + length.chrom.long[i])) {break}
          j <- j + 1
          start.pos.long[[i]][j] <- end.pos.long[[i]][j-1] + round(runif(1,tiled.lower.res,tiled.upper.res),4)
          end.pos.long[[i]][j] <- start.pos.long[[i]][j] + round(runif(1,length.clone.lower,length.clone.upper),4)
          }
          start.pos.long[[i]] <- start.pos.long[[i]][1:j-1]
          end.pos.long[[i]] <- end.pos.long[[i]][1:j-1]
  }}}


# Read in empirical distributions
# There are four empirical distributions - 
# (i) An empirical distribution of (the proportional) length of "normal"
# (2 copies) regions in the non-tiled region.

if (is.null(zerolengthnontiled))
  data(zero.length.distr.non.tiled) else {
zero.length.length.distr.non.tiled <- zerolengthnontiled}

# (ii) An empirical distribution of (the proportional) length of "normal" (2 copies) regions in the tiled region.

if (is.null(zerolengthtiled))
  data(zero.length.distr.tiled) else {
zero.length.length.distr.tiled <- zerolengthtiled}

# (iii) An empirical distribution of (the proportional) length of "non-normal" regions in the non-tiled regions.

if (is.null(nonzerolengthnontiled))
  data(non.zero.length.distr.non.tiled) else {
non.zero.length.length.distr.non.tiled <- nonzerolengthnontiled}

# (iv) An empirical distribution of (the proportional) length of "non-normal" regions in the tiled regions.

if (is.null(nonzerolengthtiled))
  data(non.zero.length.distr.tiled) else {
non.zero.length.length.distr.tiled <- nonzerolengthtiled}

# .................

# Generate mid-points based upon simulated data

  mid.point.short <- list()
  mid.point.long <- list()

  for (i in 1:22){
    mid.point.short[[i]] <- 0.5*(end.pos.short[[i]] + start.pos.short[[i]])
    mid.point.long[[i]] <- 0.5*(end.pos.long[[i]] + start.pos.long[[i]])
  }

# Function for generation of simulated aCGH data

# Generate list structure for storing simulated data

simulated.data <- list()


	
# Now generate the basic data strucutre for the chormosome. There are five strucutures for each chromosome.
# (i) the short arm

  combined.one <- combined.two <- combined.three <- combined.four <- combined.five <- list()

  for (i in 1:22){
    out.short.one <- generate.data(mid.point.short[[i]],short.tiled[i],length.short.tiled[i],end.pos.short[[i]],start.short.tiled[i], length.chrom.short[i],
                                   zero.length.distr.non.tiled = zero.length.distr.non.tiled,
                                   non.zero.length.distr.non.tiled = non.zero.length.distr.non.tiled,
                                   zero.length.distr.tiled = zero.length.distr.tiled,
                                   non.zero.length.distr.tiled = non.zero.length.distr.tiled)
    out.short.two <- generate.data(mid.point.short[[i]],short.tiled[i],length.short.tiled[i],end.pos.short[[i]],start.short.tiled[i], length.chrom.short[i],
                                   zero.length.distr.non.tiled = zero.length.distr.non.tiled,
                                   non.zero.length.distr.non.tiled = non.zero.length.distr.non.tiled,
                                   zero.length.distr.tiled = zero.length.distr.tiled,
                                   non.zero.length.distr.tiled = non.zero.length.distr.tiled)
    out.short.three <- generate.data(mid.point.short[[i]],short.tiled[i],length.short.tiled[i],end.pos.short[[i]],start.short.tiled[i], length.chrom.short[i],
                                     zero.length.distr.non.tiled = zero.length.distr.non.tiled,
                                     non.zero.length.distr.non.tiled = non.zero.length.distr.non.tiled,
                                     zero.length.distr.tiled = zero.length.distr.tiled,
                                     non.zero.length.distr.tiled = non.zero.length.distr.tiled)
    out.short.four <- generate.data(mid.point.short[[i]],short.tiled[i],length.short.tiled[i],end.pos.short[[i]],start.short.tiled[i], length.chrom.short[i],
                                    zero.length.distr.non.tiled = zero.length.distr.non.tiled,
                                    non.zero.length.distr.non.tiled = non.zero.length.distr.non.tiled,
                                    zero.length.distr.tiled = zero.length.distr.tiled,
                                    non.zero.length.distr.tiled = non.zero.length.distr.tiled)
    out.short.five <- generate.data(mid.point.short[[i]],short.tiled[i],length.short.tiled[i],end.pos.short[[i]],start.short.tiled[i], length.chrom.short[i],
                                    zero.length.distr.non.tiled = zero.length.distr.non.tiled,
                                    non.zero.length.distr.non.tiled = non.zero.length.distr.non.tiled,
                                    zero.length.distr.tiled = zero.length.distr.tiled,
                                    non.zero.length.distr.tiled = non.zero.length.distr.tiled)
    
# (ii) the long arm
    
    out.long.one <- generate.data(mid.point.long[[i]],long.tiled[i],length.long.tiled[i],end.pos.long[[i]],start.long.tiled[i], length.chrom.long[i],zero.length.distr.non.tiled = zero.length.distr.non.tiled,
         non.zero.length.distr.non.tiled = non.zero.length.distr.non.tiled,
         zero.length.distr.tiled = zero.length.distr.tiled,
         non.zero.length.distr.tiled = non.zero.length.distr.tiled)					   
    out.long.two <- generate.data(mid.point.long[[i]],long.tiled[i],length.long.tiled[i],end.pos.long[[i]],start.long.tiled[i], length.chrom.long[i],
                                  zero.length.distr.non.tiled = zero.length.distr.non.tiled,
                                  non.zero.length.distr.non.tiled = non.zero.length.distr.non.tiled,
                                  zero.length.distr.tiled = zero.length.distr.tiled,
                                  non.zero.length.distr.tiled = non.zero.length.distr.tiled)						          
    out.long.three <- generate.data(mid.point.long[[i]],long.tiled[i],length.long.tiled[i],end.pos.long[[i]],start.long.tiled[i], length.chrom.long[i],
                                    zero.length.distr.non.tiled = zero.length.distr.non.tiled,
                                    non.zero.length.distr.non.tiled = non.zero.length.distr.non.tiled,
                                    zero.length.distr.tiled = zero.length.distr.tiled,
                                    non.zero.length.distr.tiled = non.zero.length.distr.tiled)							   
    out.long.four <- generate.data(mid.point.long[[i]],long.tiled[i],length.long.tiled[i],end.pos.long[[i]],start.long.tiled[i], length.chrom.long[i],
                                   zero.length.distr.non.tiled = zero.length.distr.non.tiled,
                                   non.zero.length.distr.non.tiled = non.zero.length.distr.non.tiled,
                                   zero.length.distr.tiled = zero.length.distr.tiled,
                                   non.zero.length.distr.tiled = non.zero.length.distr.tiled)							   
    out.long.five <- generate.data(mid.point.long[[i]],long.tiled[i],length.long.tiled[i],end.pos.long[[i]],start.long.tiled[i], length.chrom.long[i],
                                   zero.length.distr.non.tiled = zero.length.distr.non.tiled,
                                   non.zero.length.distr.non.tiled = non.zero.length.distr.non.tiled,
                                   zero.length.distr.tiled = zero.length.distr.tiled,
                                   non.zero.length.distr.tiled = non.zero.length.distr.tiled)
							   
# Combining the data on the short and long arm

    combined.one[[i]] <- rbind(out.short.one,out.long.one)
    combined.two[[i]] <- rbind(out.short.two,out.long.two)
    combined.three[[i]] <- rbind(out.short.three,out.long.three)
    combined.four[[i]] <- rbind(out.short.four,out.long.four)
    combined.five[[i]] <- rbind(out.short.five, out.long.five)
  }

# Adding noise and generating the "proper" dataset

# Generate simulated data

  for (i in 1:22){
    class.matrix <- NULL
    class.matrix[[1]] <- combined.one[[i]]
    class.matrix[[2]] <- combined.two[[i]]
    class.matrix[[3]] <- combined.three[[i]]
    class.matrix[[4]] <- combined.four[[i]]
    class.matrix[[5]] <- combined.five[[i]]
	
	#true.predictions <- unlist(apply(class.matrix[[1]],1,function(x){rep(x[1],x[2])}))!=unlist(apply(class.matrix[[2]],1,function(x){rep(x[1],x[2])}))

    # Class vector
    classes <- sample(c(0,1,2,3,4),size=nArrays,replace=T,prob=c(0.2,0.2,0.2,0.2,0.2))
	
                                        # Class specific vectors of probability of each aberant segment
    class.matrix[[1]][which(class.matrix[[1]][,1]!=2),3] <- 0.7
    class.matrix[[2]][which(class.matrix[[2]][,1]!=2),3] <- 0.7
    class.matrix[[3]][which(class.matrix[[3]][,1]!=2),3] <- 0.7
    class.matrix[[4]][which(class.matrix[[4]][,1]!=2),3] <- 0.7
    class.matrix[[5]][which(class.matrix[[5]][,1]!=2),3] <- 0.7

	# Data matrix
    datamatrix <- matrix(NA,nrow=sum(class.matrix[[1]][,2]),ncol=nArrays)
    colnames(datamatrix) <- rep("Sample", nArrays)
    samples<-list()

    actual.clones1 <- unlist(apply(class.matrix[[1]],1,function(x){rep(x[1],x[2])}))
    actual.clones2 <- unlist(apply(class.matrix[[2]],1,function(x){rep(x[1],x[2])}))
    actual.clones3 <- unlist(apply(class.matrix[[3]],1,function(x){rep(x[1],x[2])}))
    actual.clones4 <- unlist(apply(class.matrix[[4]],1,function(x){rep(x[1],x[2])}))
    actual.clones5 <- unlist(apply(class.matrix[[5]],1,function(x){rep(x[1],x[2])}))
    
    class.output <- matrix(nrow=length(actual.clones1),ncol=nArrays)

	# Add individual noise and outliers to each column (sample)
    for ( t in 1:nArrays){
        # Class of current sample
      class <- classes[t]

        # Template for current class
      data <- class.matrix[[class+1]]

        # Make heterogeneous data with probabilities of aberant segments x[3] from template
		#data[,1] <- apply(data,1,function(x){ sample(c(2,x[1]),1,prob=c(1-x[3],x[3])) })
      data <- unlist(apply(data,1,function(x){rep(x[1],x[2])}))
		
        # Get proportion of normal (??tumour??)  cells (copy number 2, log2 ratio = 0)
      p <- runif(1, min=0.5, max=0.7)
		
        # Get log2 ratios given a copy number change in the 1-p proportion of tumor cells
      datamatrix[,t]<-log2((data*p+2*(1-p))/2)
      colnames(datamatrix)[t] <- paste("Sample", t)

      if (is.null(sd)) sdev <- runif(1,0.1,0.2) else {sdev <- sd}

		# Add noise 
      datamatrix[,t]<-datamatrix[,t]+rnorm(length(datamatrix[,t]),mean=0,sd=sdev)
		
		# Save parameter values
      samples[[t]] <- data.frame(p=p,log2ratios=data)

        # Write acutal classes

      class.output[,t] <- switch(class+1, actual.clones1, actual.clones2, actual.clones3, actual.clones4, actual.clones5)            
    }
    simulated.data[[i]] <- list(class.output=class.output,class.matrix=class.matrix,classes=classes,datamatrix=datamatrix,samples=samples)

  }


  for (i in 1:22){
    simulated.data[[i]]$clones <- list()
    simulated.data[[i]]$clones$start.point <- c(start.pos.short[[i]],start.pos.long[[i]])
    simulated.data[[i]]$clones$end.point <- c(end.pos.short[[i]],end.pos.long[[i]])
    simulated.data[[i]]$clones$mid.point <- c(mid.point.short[[i]],mid.point.long[[i]])
    simulated.data[[i]]$chrom <- i
  }

  SegList <- list()
  SegList$M.observed <- as.matrix(simulated.data[[1]]$datamatrix)
  SegList$M.predicted <- matrix(NA, ncol = nArrays, nrow = nrow(simulated.data[[1]]$datamatrix), dimnames = dimnames(simulated.data[[1]]$datamatrix))
  for(j in 1:nArrays){
    p = simulated.data[[1]]$samples[[j]]$p
    SegList$M.predicted[,j] <- log2((simulated.data[[1]]$samples[[j]]$log2ratios*p + 2*(1-p))/2)
  }
  SegList$state <- as.matrix(simulated.data[[1]]$class.output)
  SegList$genes <- as.data.frame(matrix(NA, nrow = length(simulated.data[[1]]$clones$mid.point), ncol = 4))
  SegList$genes[,1] <- rep(1, nrow(SegList$genes))
  SegList$genes[,2] <- simulated.data[[1]]$clones$mid.point
  SegList$genes[,3] <- simulated.data[[1]]$clones$start.point
  SegList$genes[,4] <- simulated.data[[1]]$clones$end.point
  colnames(SegList$genes) <- c("Chr", "Position", "Start", "End")
  class(SegList) = "SegList"
  
  for(i in 2:22){
    SegList.temp <- list()
    class(SegList.temp) = "SegList"
    SegList.temp$M.observed <- as.matrix(simulated.data[[i]]$datamatrix)
    SegList.temp$M.predicted <- matrix(NA, ncol = nArrays, nrow = nrow(simulated.data[[i]]$datamatrix))
    for(j in 1:nArrays){
      p = simulated.data[[i]]$samples[[j]]$p
      SegList.temp$M.predicted[,j] <- log2((simulated.data[[i]]$samples[[j]]$log2ratios*p + 2*(1-p))/2)
    }
    SegList.temp$state <- as.matrix(simulated.data[[i]]$class.output)
    SegList.temp$genes <- as.data.frame(matrix(NA, nrow = length(simulated.data[[i]]$clones$mid.point), ncol = 4))
    SegList.temp$genes[,1] <- rep(i, nrow(SegList.temp$genes))
    SegList.temp$genes[,2] <- simulated.data[[i]]$clones$mid.point
    SegList.temp$genes[,3] <- simulated.data[[i]]$clones$start.point
    SegList.temp$genes[,4] <- simulated.data[[i]]$clones$end.point
    colnames(SegList.temp$genes) <- c("Chr", "Position", "Start", "End")
  
    SegList <- rbind(SegList, SegList.temp)
  }
  SegList$design <- rep(1, nArrays)
  SegList
}    


# Function for generating the points on the short arm.

"generate.data" <-
function(mid.point,
         tiled,
         length.tiled,
         end.pos,
         start.tiled,
         length.chrom,
         zero.length.distr.non.tiled. = zero.length.distr.non.tiled,
         non.zero.length.distr.non.tiled. = non.zero.length.distr.non.tiled,
         zero.length.distr.tiled. = zero.length.distr.tiled,
         non.zero.length.distr.tiled. = non.zero.length.distr.tiled) {
  data<-numeric()
  if (tiled == 0){
    j <- 0
  repeat {
    if (j >= length(mid.point)) {break}
	{
 	  state <- sample(c(0,1,2,3,4,5),1,prob=c(0.04,0.15,0.5,0.15,0.10,0.06))
 	  if (state==2){ l <-
  max(1,floor(sample(zero.length.distr.non.tiled. , 1))*length(mid.point))}
 	  else {l <-
  max(1,floor(sample(non.zero.length.distr.non.tiled. , 1)*length(mid.point)))
  }
 	  if (j+l > length(mid.point)) { l <- length(mid.point) -
  j }
  	  data<-rbind(data,c(state,l,1))
  	  j<-j+l
	}}} else {
      if (tiled == 1){
        if (length.tiled != length.chrom){
          j <- 0
        repeat {
          if (j >= length(mid.point[mid.point < start.tiled])) {break}
          {
            state <- sample(c(0,1,2,3,4,5),1,prob=c(0.04,0.15,0.5,0.15,0.10,0.06))
      	    if (state==2){ l <-
                          max(floor(sample(zero.length.distr.non.tiled. , 1)*(length(mid.point[mid.point < start.tiled]))),1)} else
            {l <-
  max(1,floor(sample(non.zero.length.distr.non.tiled. , 1)*(length(mid.point[mid.point
  < start.tiled]))))}
            if (j+l > (length(mid.point[mid.point <
  start.tiled]))) { l <- length(mid.point[mid.point <
  start.tiled]) - j}
  	        data<-rbind(data,c(state,l,1))
 	        j<-j+l}}
            j <- length(mid.point[mid.point < start.tiled])
         repeat {
            if (j >= length(mid.point[mid.point < (start.tiled + length.tiled)])) {break}
            {
            state <- sample(c(0,1,2,3,4,5),1,prob=c(0.04,0.15,0.5,0.15,0.10,0.06))
      	    if (state==2){ l <-
                          max(1,floor(sample(zero.length.distr.tiled. , 1)*(length(mid.point[mid.point
  > start.tiled & mid.point < (start.tiled + length.tiled)]))))} else
            {l <-
  max(1,floor(sample(non.zero.length.distr.tiled. , 1)*(length(mid.point[mid.point >
  start.tiled & mid.point < (start.tiled +
  length.tiled)]))))}
            if (j+l > (length(mid.point[mid.point < (start.tiled +
  length.tiled)]))) { l <- length(mid.point[mid.point <
  (start.tiled + length.tiled)]) - j}
  	        data<-rbind(data,c(state,l,1))
 	        j<-j+l}}
          j <- length(mid.point[mid.point < (start.tiled+ length.tiled)])
          repeat {
            if (j >= length(mid.point)) {break}
            {
             state <- sample(c(0,1,2,3,4,5),1,prob=c(0.04,0.15,0.5,0.15,0.10,0.06))
      	    if (state==2){ l <-
                          max(1,floor(sample(zero.length.distr.non.tiled. , 1)*(length(mid.point[mid.point
  >  (start.tiled + length.tiled)]))))} else
            {l <-
  max(1,floor(sample(non.zero.length.distr.non.tiled. , 1)*(length(mid.point[mid.point
  > (start.tiled + length.tiled)]))))}
            if (j+l > (length(mid.point))) { l <-
  length(mid.point) - j}
  	        data<-rbind(data,c(state,l,1))
 	        j<-j+l}}} else {
              if (length.tiled == length.chrom){
              j <- 0
              repeat {
          if (j >= length(mid.point)) {break}
          {
            state <- sample(c(0,1,2,3,4,5),1,prob=c(0.04,0.15,0.5,0.15,0.10,0.06))
            if (state==2){ l <- max(1,floor(sample(zero.length.distr.tiled. , 1)*length(mid.point)))}
 	        else {l <-
  max(1,floor(sample(non.zero.length.distr.tiled. , 1)*length(mid.point)))
  }
 	        if (j+l > length(mid.point)) { l <-
  length(mid.point) - j}
  	        data<-rbind(data,c(state,l,1))
 	        j<-j+l
	}}}}}}
	return(data)}
