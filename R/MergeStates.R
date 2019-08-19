"mergeStates" <- function(segList, MergeType = 1, pv.thres=0.0001, ansari.sign=0.01, minDiff = 0.25){
  if (MergeType == 1) {
    MergeLevels.new(segList, pv.thres = pv.thres, ansari.sign = ansari.sign)
  }
  else if (MergeType == 2) {
    MergeLevels.old(segList, minDiff = minDiff)
  }
  else { 
    stop("\nMergeType argument must be either 1 or 2.  \nPlease see the help file for more details")
  }
}

########################


# Level Merging Function

MergeLevels.new <- function(seg.info,pv.thres=0.0001,ansari.sign=0.01){
  chrom <- seg.info$genes$Chr
  for(k in 1:length(seg.info$M.observed[1,])) {
    counter = 0
    for(n in 1:length(unique(chrom))) {

      vecObs <- seg.info$M.observed[chrom == n,k]
      vecPred <- seg.info$M.predicted[chrom == n, k]
      states <- seg.info$state[chrom == n,k]
      num.states <- length(unique(states))
      if (num.states > 1) {
        #threshold for wilcoxon test
        pv.thres=0.0001


        sq<-numeric()  # Initializing threshold vector for keeping track of thresholds
        j=0   #initializing threshold index (threshold count)

                                        #initializing ansari p-values
        ansari=numeric()
        ansari.sign=0.05

# Initialize levels count
        lv=rep(NA, length(sq))

# Start with threshold 0.01, and flag=0 indicating significance not yet reached, backtrack not begun
        thresAbs=0.05
        flag=0
        while (1){
          j=j+1
                                        # Save threshold
          sq[j]<-thresAbs

	# temporary predicted values (to be updated)
          vecPredNow=vecPred
	#unmissing unique segment medians
          mnNow=unique(vecPred)
          mnNow=mnNow[!is.na(mnNow)]

	#continuing indicator otherwise get out of the loop
          cont=0

          while(cont==0 & length(mnNow)>1) {

            mnNow=sort(mnNow) #currennt sorted vector of means
       	#dst=as.matrix(dist(mnNow)) #distances between all pairs of segment means
            dst=as.matrix(dist(2*2^mnNow))

		#To extract indices of the mean pairs:
            dst.row=matrix(1:nrow(dst), nrow=nrow(dst), ncol=ncol(dst))
            dst.col=matrix(1:ncol(dst), nrow=nrow(dst), ncol=ncol(dst), byrow=TRUE)
            lt=lower.tri(dst)  #lower triangle of the distance matrix
            dst=matrix(cbind(dst[lt], dst.row[lt], dst.col[lt]),ncol=3)

        # Remove any pairs more than 1 level apart
            dst=matrix(dst[abs(dst[,2]-dst[,3])<=1,],ncol=3)

		# Scale differences in mean by max median of the two observed vectors
		#dst[,1]=dst[,1]*apply(cbind(abs(mnNow[dst[,2]]),abs(mnNow[dst[,3]])),1,max)

		#order distance between means with the closest on top and corresponding indices
            dst=matrix(dst[order(dst[,1]),],ncol=3)

		#for each pair of means
            for (i in 1:nrow(dst)) 	{
			#set continuity index to "NOT continue" (=1)
              cont=1
			#test for combining of the two segment means
              out=combine.func(diff=dst[i,1],vecObs, vecPredNow, mnNow, mn1=mnNow[dst[i,2]], mn2=mnNow[dst[i,3]], pv.thres=pv.thres, thresAbs=thresAbs)
			#if combine?
              if (out$pv > pv.thres) {

               	#set continuity index to "YES" (=0) and break out of the current pairs loop
                cont=0

				#update predicted values and segments
                vecPredNow=out$vecPredNow
                mnNow=out$mnNow
                break
              }		
            }		
          }
	# When done merging for a given threshold, test for significance
          ansari[j]=ansari.test(sort(vecObs-vecPredNow), sort(vecObs-vecPred))$p.value
          lv[j]=length(mnNow) # get number of levels

	# If p.value is less than the significance threshold, set backtracking flag=1 (backtracking on)
          if(ansari[j]<ansari.sign){
            flag=1
          }
          if(flag==2){ break } 
	
	# If backtracking is on, a smaller threshold is attempted
          if (flag){

		# If backtracking is on and p.value is higher than sign threshold, stop
            if (ansari[j]>ansari.sign | thresAbs == 0){
            # Don't merge at all if all tested threshold including 0 is significant
              if (ansari[j] <= ansari.sign) {
                vecPredNow=vecPred
                mnNow=unique(vecPred)
                mnNow=mnNow[!is.na(mnNow)]
              }
              break
            }  

        # Attempt smaller threshold
            else {thresAbs=signif(thresAbs-0.005,3) }	
          }
          else {thresAbs=thresAbs+0.1} # Increase threshold if backtracking is not on

	# Control so it won't keep going, max threshold = 1 and if sign not reached, threshold = 0
          if (thresAbs >= 1){
            thresAbs=0
            flag=2
          }
        }

# Return list of results

        seg.info$M.predicted[chrom == n, k] <- vecPredNow
        seg.info$num.states[n,k] <- length(unique(vecPredNow))

        states.uniq <- unique(vecPredNow)
        temp.states <- c(1:length(vecPredNow))

        for(g in 1:length(unique(vecPredNow))) {
          temp.states[vecPredNow == states.uniq[g]] <- g
        }
        seg.info$state[chrom == n, k] <- temp.states
      }     
    }
  }
  new("SegList",seg.info)
}

#This is the original merge function that is supplied with the aCGH library

"MergeLevels.old" <- 
  function(segList, minDiff = 0.25){
    chrom <- segList$genes$Chr
    for(i in 1:ncol(segList)){
      for(j in 1:length(unique(chrom))){
        
        ind.nonna <- which(!is.na(segList$M.observed[ ,i]))
        segList <- segList[ind.nonna,]
        
        statesr <- list(segList$state[chrom == j,i],
			segList$M.predicted[chrom == j,i],
			segList$dispersion[chrom == j,i],
			segList$M.observed[chrom == j,i])
		#	segList[[5]][chrom == j,i],
		#	segList[[6]][chrom == j,i])
			#probably should check and remove NA's here.
                                        #need to check and see if there ever could be any
        states <- statesr[[1]]
        M.predicted <- statesr[[2]]
        M.observed <- statesr[[4]]
        num.states <- length(unique(states))
        if (num.states > 1) {
          for (m in 1:(num.states - 1)) {
            states.uniq <- unique(states)
            pred.states.uniq <- rep(0, length(states.uniq))
            for(s in 1:length(states.uniq)) {
              M.predicted[states == states.uniq[s]] <- median(M.observed[states == states.uniq[s]])
              pred.states.uniq[s] <- (M.predicted[states == states.uniq[s]])[1]
            }
            dst <- abs(dist(pred.states.uniq))
            if (min(dst) >= minDiff){
              segList$state[chrom == j, i] <- states
              segList$M.predicted[chrom ==j, i] <- M.predicted
              break
            }
            else {
              pred.dist.matr <- as.matrix(dst)
              for (s1 in 1:(nrow(pred.dist.matr) - 1)) {
                for (s2 in (s1 + 1):ncol(pred.dist.matr)) {
                  if (pred.dist.matr[s1,s2] == min(dst)) {
                    states[states == states.uniq[s2]] <- states.uniq[s1]
                    M.predicted[states == states.uniq[s1]] <- median(M.observed[states == states.uniq[s1]])
                    break
                  }
                }
              }
            }
          }
        }
        segList$state[chrom == j, i] <- states
        segList$M.predicted[chrom ==j, i] <- M.predicted
        segList$num.states[j,i] <- length(unique(states))
      }
    }
    segList
  }



