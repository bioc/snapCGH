"find.genomic.events" <-
function (segList, maxChrom = 22, factor = 5, maxClones = 1, 
    maxLen = 1000, absValSingle = 1, absValRegion = 1, diffVal1 = 1, 
    diffVal2 = 0.5, maxSize = 10000, pChrom.min = 0.9, medChrom.min = 0.1,
    #The following are used in the computeSD.samples function
          maxmadUse = 0.2, maxmedUse = 0.2, 
    maxState = 3, maxStateChange = 10, minClone = 5) 
{
 #computing Std. Devs of aCGH samples
 sd.samples <- computeSD.func(segList, maxmadUse, maxmedUse, maxState, maxStateChange, minClone, maxChrom)
 
    l2r <- log2.ratios(segList)
    genes <- segList$genes
    
    ncols <- ncol(l2r)
    cat("Finding outliers\n")
    outliers <- findOutliers.func(segList, thres = sd.samples$madGenome, 
        factor = factor)
    cat("Finding focal low level aberrations\n")
    aberrations <- findAber.func(segList, maxClones = maxClones, maxLen = maxLen)
    cat("Finding transitions\n")
    transitions <- findTrans.func(segList, outliers = outliers$outlier, 
        aber = aberrations$aber)
    cat("Finding focal amplifications\n")
    amplifications <- findAmplif.func(segList, absValSingle = absValSingle, 
        absValRegion = absValRegion, diffVal1 = diffVal1, diffVal2 = diffVal2, 
        maxSize = maxSize, translen.matr = transitions$translen.matrix, 
        trans.matr = transitions$trans.matr, aber = aberrations$aber, 
        outliers = outliers$outlier, pred = outliers$pred.out, 
        pred.obs = outliers$pred.obs.out)
    numTrans <- matrix(0, nrow = maxChrom, ncol = ncols)
    numAmplif <- matrix(0, nrow = maxChrom, ncol = ncols)
    numAber <- matrix(0, nrow = maxChrom, ncol = ncols)
    numOutlier <- matrix(0, nrow = maxChrom, ncol = ncols)
    numTrans.binary <- matrix(0, nrow = maxChrom, ncol = ncols)
    numAmplif.binary <- matrix(0, nrow = maxChrom, ncol = ncols)
    numAber.binary <- matrix(0, nrow = maxChrom, ncol = ncols)
    numOutlier.binary <- matrix(0, nrow = maxChrom, ncol = ncols)
    wholeChromGainLoss <- matrix(0, nrow = maxChrom, ncol = ncols)
    sizeAmplicon <- numAmplicon <- matrix(0, nrow = maxChrom, 
        ncol = ncols)
    propClones <- matrix(0, nrow = maxChrom, ncol = ncols)
    pvClones <- matrix(0, nrow = maxChrom, ncol = ncols)
    medClones <- matrix(0, nrow = maxChrom, ncol = ncols)
    p.min <- pChrom.min
    pv.max <- 1e-04
    med.min <- medChrom.min
    chr <- segList$genes$Chr
    kb <- segList$genes$Position
    for (j in 1:maxChrom) {
        cat("Processing chromosome ", j, "\n")
        ind <- chr == j
        trans <- transitions$trans.matrix[ind, , drop = FALSE]
        amplif <- amplifications$amplif[ind, , drop = FALSE]
        aber <- aberrations$aber[ind, , drop = FALSE]
        outlier <- outliers$outlier[ind, , drop = FALSE]
        for (i in 1:ncols) {
            numTrans[j, i] <- sum(trans[, i] == 1, na.rm = TRUE)
            if (numTrans[j, i] > 0) 
                numTrans.binary[j, i] <- 1
            else {
                obs <- l2r[ind, i]
                obs <- obs[aber[, i] == 0 & outlier[, i] == 0]
                obs <- na.omit(obs)
                if (length(obs) != 0) {
                  p <- if (median(obs) >= 0) 
                    mean(obs >= 0)
                  else mean(obs < 0)
                  propClones[j, i] <- p
                  pv <- 1 - pnorm((p - 0.5)/sqrt((0.5^2)/length(obs)))
                  pvClones[j, i] <- pv
                  medClones[j, i] <- median(obs)
                  if ((p >= p.min) && (abs(median(obs)) >= med.min)) {
                    if (pv <= pv.max) 
                      wholeChromGainLoss[j, i] <- if (median(obs) >= 
                        0) 
                        1
                      else -1
                  }
                  else wholeChromGainLoss[j, i] <- 0
                }
                else wholeChromGainLoss[j, i] <- 0
            }
            numAmplif[j, i] <- sum(amplif[, i] == 1, na.rm = TRUE)
            if (numAmplif[j, i] > 0) 
                numAmplif.binary[j, i] <- 1
            numAber[j, i] <- sum(aber[, i] == 1, na.rm = TRUE)
            if (numAber[j, i] > 0) 
                numAber.binary[j, i] <- 1
            numOutlier[j, i] <- sum(outlier[, i] == 1, na.rm = TRUE)
            if (numOutlier[j, i] > 0) 
                numOutlier.binary[j, i] <- 1
            amp <- amplif[, i]
            ind.na <- which(is.na(amp))
            amp <- amp[-ind.na]
            try1 <- diff(amp)
            tmps <- which(try1 == 1) + 1
            tmpe <- which(try1 == -1)
            if (length(tmps) > length(tmpe)) 
                tmpe <- c(tmpe, length(ind))
            if (length(tmps) < length(tmpe)) 
                tmps <- c(1, tmps)
            if (length(tmpe) == length(tmps)) {
                kb.ind <- kb[ind][-ind.na]
                tmplen <- (kb.ind[tmpe] - kb.ind[tmps]) + rep(1000, 
                  length(tmpe))
                sizeAmplicon[j, i] <- sum(tmplen)
                numAmplicon[j, i] <- length(tmpe)
            }
        }
    }
    ge <- list(num.transitions = numTrans, num.amplifications = numAmplif, 
        num.aberrations = numAber, num.outliers = numOutlier, 
        num.transitions.binary = numTrans.binary, num.amplifications.binary = numAmplif.binary, 
        num.aberrations.binary = numAber.binary, num.outliers.binary = numOutlier.binary, 
        whole.chrom.gain.loss = wholeChromGainLoss, size.amplicons = sizeAmplicon, 
        num.amplicons = numAmplicon, outliers = outliers, aberrations = aberrations, 
        transitions = transitions, amplifications = amplifications)
    new("GEList", ge)    
}

"findAmplif.func" <-
function (segList, absValSingle = 1, absValRegion = 1.5, diffVal1 = 1, 
    diffVal2 = 0.5, maxSize = 10000, translen.matr = res3$translen.matrix, 
    trans.matr = res3$trans.matr, aber = res2$aber, outliers = res1$outlier, 
    pred = res1$pred.out, pred.obs = res1$pred.obs.out) 
{
    chrom <- segList$genes$Chr
    kb <- segList$genes$Position
    amplif.matrix <- matrix(0, nrow = length(kb), ncol = ncol(segList))
    for (i in 1:ncol(segList)) {
        for (j in 1:length(unique(chrom))) {
            ind.nonna <- (1:length(segList$M.observed[chrom == j, i]))[!is.na(segList$M.observed[chrom == 
                j, i])]
            abernow <- aber[chrom == j, i][ind.nonna]
            outliersnow <- outliers[chrom == j, i][ind.nonna]
            prednow <- pred[chrom == j, i][ind.nonna]
            predobsnow <- pred.obs[chrom == j, i][ind.nonna]
            obsnow <- segList$M.observed[chrom == j, i][ind.nonna]
            transnow <- trans.matr[chrom == j, i][ind.nonna]
            translennow <- translen.matr[chrom == j, i][ind.nonna]
            amplifnow <- rep(0, length(obsnow))
            amplifnow[outliersnow == 1 & ((obsnow - prednow) >= 
                diffVal1)] <- 1
            amplifnow[outliersnow == 1 & ((obsnow - prednow) >= 
                diffVal2) & obsnow >= absValSingle] <- 1
            indaber <- (1:length(amplifnow))[abernow == 1]
            if (length(indaber) > 0) {
                indstretch <- (1:length(amplifnow))[abernow == 
                  0 & outliersnow == 0]
                for (m in 1:length(indaber)) {
                  stretchleft <- max(0, max(indstretch[indstretch < 
                    indaber[m]]), na.rm = TRUE)
                  stretchright <- min((length(amplifnow) + 1), 
                    min(indstretch[indstretch > indaber[m]]), 
                    na.rm = TRUE)
                  if (stretchleft == 0) {
                    mx <- prednow[stretchright]
                  }
                  else if (stretchright == (length(amplifnow) + 
                    1)) {
                    mx <- prednow[stretchleft]
                  }
                  else {
                    mx <- max(prednow[stretchleft], prednow[stretchright])
                  }
                  if (!is.na(mx)) {
                    if (((predobsnow[indaber[m]] - mx) >= diffVal1) || 
                      ((predobsnow[indaber[m]] - mx) >= diffVal2 && 
                        (predobsnow[indaber[m]] >= absValSingle))) {
                      amplifnow[indaber[m]] <- 1
                    }
                  }
                }
            }
            amplifnow[abernow == 0 & outliersnow == 0 & obsnow >= 
                absValRegion & translennow <= maxSize] <- 1
            amplif.matrix[chrom == j, i][ind.nonna] <- amplifnow
            amplif.matrix[chrom == j, i][-ind.nonna] <- NA
        }
    }
    list(amplif = amplif.matrix)
}

"findTrans.func" <-
function (segList, outliers = res1$outliers, aber = res2$aber) 
{
    chrom <- segList$genes$Chr
    kb <- segList$genes$Position
    trans.matrix <- matrix(0, nrow = length(chrom), ncol = ncol(segList))
    translen.matrix <- matrix(0, nrow = length(chrom), ncol = ncol(segList))
    for (i in 1:ncol(segList)) {
        for (j in 1:length(unique(chrom))) {
            ind.nonna <- (1:length(segList$M.observed[chrom == j, i]))[!is.na(segList$M.observed[chrom == 
                j, i])]
            kbnow <- kb[chrom == j][ind.nonna]
            states <- segList$state[chrom == j, i][ind.nonna]
            outliersnow <- outliers[chrom == j, i][ind.nonna]
            abernow <- aber[chrom == j, i][ind.nonna]
            transnow <- rep(0, length(states))
            translennow <- rep(0, length(states))
            states.diff <- diff(states[abernow == 0])
            ind <- (1:length(states.diff))[states.diff != 0]
            if (length(ind) > 0) {
                start <- ind + 1
                end <- ind
                transnow[abernow == 0][start] <- 1
                transnow[abernow == 0][end] <- 2
            }
            st <- c(1, (1:length(transnow))[transnow == 1])
            en <- c((1:length(transnow))[transnow == 2], length(transnow))
            for (m in 1:length(st)) {
                translennow[st[m]:en[m]] <- kbnow[en[m]] - kbnow[st[m]]
            }
            translen.matrix[chrom == j, i][ind.nonna] <- translennow
            transnow[abernow == 1] <- 3
            trans.matrix[chrom == j, i][ind.nonna] <- transnow
            trans.matrix[chrom == j, i][-ind.nonna] <- NA
            translen.matrix[chrom == j, i][-ind.nonna] <- NA
        }
    }
    list(trans.matrix = trans.matrix, translen.matrix = translen.matrix)
}
"findAber.func" <-
function (segList, maxClones = 1, maxLen = 1000) 
{
    chrom <- segList$genes$Chr
    kb <- segList$genes$Position
    aber <- matrix(0, nrow = length(chrom), ncol = ncol(segList))
    for (i in 1:ncol(segList)) {
        for (j in 1:length(unique(chrom))) {
            st1 <- segList$M.observed[chrom == j, i]
            ind.nonna <- which(!is.na(st1))
            states <- segList$state[chrom == j, i][ind.nonna]
            kbnow <- kb[chrom == j][ind.nonna]
            abernow <- rep(0, length(kbnow))
            num <- 1
            for (m in 2:length(states)) {
                if (states[m - 1] != states[m]) {
                  if (m == 2) {
                    abernow[1] <- 1
                  }
                  if (m == 3) {
                    abernow[1:2] <- 1
                  }
                  if (m == length(states)) {
                    abernow[length(states)] <- 1
                  }
                  if (m == (length(states) - 1)) {
                    abernow[(length(states) - 1):(length(states))] <- 1
                  }
                  if (m <= length(states)) {
                    if ((num <= maxClones) || ((kbnow[m - 1] - 
                      kbnow[m - num]) <= maxLen)) {
                      abernow[(m - num):(m - 1)] <- 1
                    }
                  }
                  num <- 1
                }
                else {
                  num <- num + 1
                }
            }
            aber[chrom == j, i][ind.nonna] <- abernow
            aber[chrom == j, i][-ind.nonna] <- NA
        }
    }
    list(aber = aber)
}

"findOutliers.func" <-
function (segList, thres = madGenome, factor = 4) 
{
    thres <- thres * factor
    chrom <- segList$genes$Chr
    outlier <- matrix(0, nrow = length(chrom), ncol = ncol(segList))
    states.out <- matrix(0, nrow = length(chrom), ncol = ncol(segList))
    pred.out <- matrix(0, nrow = length(chrom), ncol = ncol(segList))
    pred.obs.out <- matrix(0, nrow = length(chrom), ncol = ncol(segList))
    for (i in 1:ncol(segList)) {
        for (j in 1:length(unique(chrom))) {
            ind.nonna <- which(!is.na(segList$M.observed[chrom == j,i]))
            obs <- segList$M.observed[chrom == j, i][ind.nonna]
            pred <- segList$M.predicted[chrom == j, i][ind.nonna]
            pred.obs <- pred
#            states <- matrix(0, nrow = length(obs), ncol = 1)
            states <- segList$state[chrom ==j, i][ind.nonna]

            for (k in 1:length(obs))
              {
                md <-  median(obs[states == states[k]],na.rm = TRUE)
                if ((obs[k] >  md + thres[i]) || (obs[k] < md - thres[i]))
                {
                    outlier[chrom==j, i][ind.nonna][k] <- 1
                    ##assigning observed value for a clone instead of predicted
                    pred.obs[k] <- obs[k]
                }
            }
            states.uniq <- unique(states)
            for (m in 1:length(states.uniq))
            {
                
                ##predictions for all
                pred[states==states.uniq[m]] <-
                    median(obs[states == states.uniq[m] & outlier[chrom==j, i][ind.nonna] == 0])
                ##predictions for non-outliers only
                pred.obs[states == states.uniq[m] & outlier[chrom==j, i][ind.nonna] == 0] <-
                    median(obs[states==states.uniq[m] & outlier[chrom==j, i][ind.nonna] == 0])
                
              }            
            pred.obs.out[chrom==j, i][ind.nonna] <- pred.obs
            pred.out[chrom==j, i][ind.nonna] <- pred
            outlier[chrom==j, i][-ind.nonna] <- NA
            pred.obs.out[chrom==j, i][-ind.nonna] <- NA
            pred.out[chrom==j, i][-ind.nonna] <- NA
            
        }
    }
    list(outlier=outlier, pred.obs.out=pred.obs.out, pred.out=pred.out)
}

"computeSD.func" <-
function (segList, maxmadUse = 0.2, maxmedUse = 0.2, 
    maxState = 3, maxStateChange = 10, minClone = 5, maxChrom = 22) 
{
    chrom <- segList$genes$Chr
    madChrom <- matrix(NA, nrow = length(unique(chrom)), ncol = length(segList$M.predicted[1,]))
    madGenome <- rep(NA, length(segList$M.predicted[1,]))
    for (i in 1:length(segList$M.predicted[1,])) {
        mad.tmp <- NA
        for (j in 1:maxChrom) {
        
        # this is some checking for NA's in the data
        #can there really be any since we've removed some and imputed the rest in earlier stages????
        
            ind.nonna <- (1:length(segList$M.observed[chrom == j, i]))[!is.na(segList$M.observed[chrom == 
               j, i])]
            
            mad.tmp1 <- NA
            states.uniq <- unique(segList$state[chrom == j, i])
            states <- segList$state[chrom == j, i]
            obs <- segList$M.observed[chrom == j, i]
            states.change <- diff(states)
            states.change.num <- length(states.change[states.change != 0])
            if ((length(states.uniq) <= maxState) && (states.change.num <= maxStateChange)) {
                for (k in states.uniq) {
                  obs.state <- obs[states == k]
                  md <- mad(obs.state, na.rm = TRUE)
                  med <- median(obs.state, na.rm = TRUE)
                  if ((length(obs.state) >= minClone) && (md <= 
                    maxmadUse) && (abs(med) <= maxmedUse)) {
                    mad.tmp1 <- c(mad.tmp1, md)
                  }
                }
            }
            if (length(mad.tmp1[!is.na(mad.tmp1)]) > 0) {
                mad.tmp1 <- mad.tmp1[-1]
                mad.tmp <- c(mad.tmp, mad.tmp1)
                madChrom[j, i] <- median(mad.tmp1)
            }
        }
        if (length(mad.tmp[!is.na(mad.tmp)]) > 0) {
            mad.tmp <- mad.tmp[-1]
            madGenome[i] <- median(mad.tmp)
        }
    }
    if (length(madGenome[is.na(madGenome)]) > 0) 
        cat("Warning: MAD may not have been computed for one of the\nsamples\n")
    list(madChrom = madChrom, madGenome = madGenome)
}
