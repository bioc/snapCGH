fractionAltered <- function(input, thres=0.25, factor=2.5, samplenames = colnames(input$M), chrominfo=chrominfo.basepair) {

#        #check if sd.samples are non-empty:
	#if (!is.null(sd.samples(aCGH.obj)) && (factor > 0)) {
#		thres <- factor*(sd.samples(aCGH.obj)$madGenome)
#	}
    
	data <- log2ratios(input)
	data.thres <- threshold(data, thresAbs=thres)
	datainfo <- input$genes
	ord <- order(datainfo$Chr,datainfo$Position)
	datainfo <- datainfo[ord,]

	nClones <- nrow(datainfo)
	prevChrom <- 0
	cloneLens <- rep(0,nClones)
	for (i in 1:nClones) {
		nextChrom <- datainfo$Chr[i+1]
		if (i==nClones) nextChrom <- 9999
		startPos <- ifelse(datainfo$Chr[i]==prevChrom,sum(datainfo$Position[i-1],datainfo$Position[i])/2,0)
		endPos <- ifelse(datainfo$Chr[i]==nextChrom,
					sum(datainfo$Position[i],datainfo$Position[i+1])/2,
					chrominfo$length[datainfo$Chr[i]])
		cloneLens[i] <- endPos-startPos
		prevChrom <- datainfo$Chr[i]
	}
	cloneLens <- cloneLens[order(ord)]
	loss <- gain <- rep(0, ncol(input))
	for (i in 1:ncol(input)) {
		clones.ind.na <- which(!is.na(data.thres[,i]))
		totCloneLens <- sum(cloneLens[clones.ind.na], na.rm = TRUE)
		gain[i] <- sum(cloneLens[data.thres[clones.ind.na,i]==1], na.rm = TRUE)/totCloneLens #proportion gained
		loss[i] <- sum(cloneLens[data.thres[clones.ind.na,i]==-1], na.rm = TRUE)/totCloneLens #proportion lost
	}
	list(gainP = gain, lossP = loss)
}

"plotGainLoss" <-
function (segList, resT = NULL, pheno = rep(1, ncol(log2ratios(segList))), 
    chrominfo = chrominfo.basepair, X = TRUE, Y = FALSE, 
    rsp.uniq = unique(pheno), all = length(rsp.uniq) == 1 && 
        is.null(resT), titles = if (all) "All Samples" else rsp.uniq, 
    cutplot = 0, thres = 0.25, factor = 2.5, ylm = c(-1, 1), 
    p.thres = c(0.01, 0.05, 0.1), numaut = 22, onepage = TRUE, 
    colored = TRUE) 
{

  sd.samples <- computeSD.func(segList)
  thres <- factor * sd.samples$madGenome
  
    col.scheme <- if (colored) 
        list(pal = c("red", "blue", "green", "orange")[1:length(p.thres)], 
            gain.low = "white", gain.high = "green", loss.low = "red", 
            loss.high = "white", abline1 = "blue", abline2 = "grey50", 
            mtext = "red", kb.loc = "blue", abline3 = "black", 
            abline4 = "grey50", )
    else list(pal = c("grey10", "grey40", "grey70", "grey90")[1:length(p.thres)], 
        gain.low = "grey50", gain.high = "grey0", loss.low = "grey0", 
        loss.high = "grey50", abline1 = "grey50", abline2 = "grey50", 
        mtext = "black", kb.loc = "black", abline3 = "black", 
        abline4 = "grey50", )
    data <- log2ratios(segList)
    datainfo <- segList$genes
    rsp.uniq <- sort(rsp.uniq)
    colmatr <- if (length(rsp.uniq) > 1) 
        t(sapply(rsp.uniq, function(rsp.uniq.level) ifelse(pheno == 
            rsp.uniq.level, 1, 0)))
    else matrix(rep(1, length(pheno)), ncol = length(pheno), 
        nrow = 1)
    nr <- nrow(colmatr)
    if (!is.null(resT)) 
        nr <- nr + 1
    tmp <- as.data.frame(matrix(0, ncol = 2, nrow = 1))
    colnames(tmp) <- c("gainP", "lossP")
    gainloss <- lapply(1:nrow(colmatr),
                       function(j)
                       gainLoss(input = data,
                                cols = which(colmatr[j, ] == 1),
                                thres = thres)
                       )  
    dt <- data[, colmatr[1, ] == 1, drop = FALSE]
    rsp <- rep(1, ncol(dt))
    if (nrow(colmatr) > 1) 
        for (j in 2:nrow(colmatr)) {
            dt <- cbind(dt, data[, colmatr[j, ] == 1])
            rsp <- c(rsp, rep(j, sum(colmatr[j, ] == 1)))
        }
    rsp <- rsp - 1
    numchr <- numaut
    if (X) 
        numchr <- numchr + 1
    if (Y) 
        numchr <- numchr + 1
    chrominfo <- chrominfo[1:numchr, ]
    start <- c(0, cumsum(chrominfo$length))
    kb.loc <- datainfo$Position
    for (i in 1:nrow(chrominfo)) kb.loc[datainfo$Chr == i] <- start[i] + 
        datainfo$Position[datainfo$Chr == i]
    chrom.start <- c(0, cumsum(chrominfo$length))[1:numchr]
    chrom.centr <- chrom.start + chrominfo$centr
    chrom.mid <- chrom.start + chrominfo$length/2
    par(mfrow = c((if (onepage) nr else 1), 1), lab = c(1, 8, 
        7), tcl = -0.2, xaxs = "i")
    for (g in 1:length(titles)) {
        gl <- gainloss[[g]]
        tl <- as.character(titles[g])
        ylm[1] <- min(ylm, min(gl$lossP))
        ylm[2] <- max(ylm, max(gl$gainP))
        ind <- which(gl$gainP >= cutplot)
        if (colored) {
            plot(kb.loc[ind], gl$gainP[ind], col = "green", type = "h", 
                xlab = "chromosome number", ylab = "Fraction gained or lost", 
                pch = 18, main = tl, ylim = ylm, xlim = c(0, 
                  max(cumsum(chrominfo$length), kb.loc[ind], 
                    rm.na = TRUE)), xaxt = "n")
            axis(side = 1, at = kb.loc[ind][1], label = "", tick = FALSE)
            ind <- gl$lossP >= cutplot
            points(kb.loc[ind], -gl$lossP[ind], col = "red", 
                type = "h")
        }
        else {
            plot(kb.loc[ind], gl$gainP[ind], col = "grey10", 
                type = "h", xlab = "chromosome number", ylab = "Fraction gained or lost", 
                pch = 18, main = tl, ylim = ylm, xlim = c(0, 
                  max(cumsum(chrominfo$length), kb.loc[ind], 
                    rm.na = TRUE)), xaxt = "n")
            axis(side = 1, at = kb.loc[ind][1], label = "", tick = FALSE)
            ind <- gl$lossP >= cutplot
            points(kb.loc[ind], -gl$lossP[ind], col = "grey50", 
                type = "h")
        }
        abline(h = 0)
        abline(v = cumsum(chrominfo$length), col = col.scheme$abline1)
        abline(v = chrom.centr, lty = 2, col = col.scheme$abline2)
        for (i in seq(2, numaut, b = 2)) mtext(paste("", i), 
            side = 3, at = (chrom.mid[i]), line = 0.3, col = col.scheme$mtext, 
            cex.main = 0.5)
        for (i in seq(1, numaut, b = 2)) mtext(paste("", i), 
            side = 1, at = (chrom.mid[i]), line = 0.3, col = col.scheme$mtext, 
            cex.main = 0.5)
        if (X) 
            if (is.even(numaut)) 
                mtext("X", side = 1, at = (chrom.mid[numaut + 
                  1]), line = 0.3, col = col.scheme$mtext, cex.main = 0.5)
            else mtext("X", side = 3, at = (chrom.mid[numaut + 
                1]), line = 0.3, col = col.scheme$mtext, cex.main = 0.5)
        if (Y) 
            if (is.even(numaut)) 
                mtext("Y", side = 3, at = (chrom.mid[numaut + 
                  2]), line = 0.3, col = col.scheme$mtext, cex.main = 0.5)
            else mtext("Y", side = 1, at = (chrom.mid[numaut + 
                2]), line = 0.3, col = col.scheme$mtext, cex.main = 0.5)
    }
    if (!is.null(resT)) {
        res <- resT[order(resT$index), ]
        maxT <- res$adjp
        teststat <- abs(res$teststat)
        st <- sapply(p.thres, function(threshold) {
            if (any(maxT <= threshold)) 
                min(teststat[maxT <= threshold])
            else NA
        })
        plot(kb.loc, teststat, col = col.scheme$kb.loc, ylim = c(0, 
            max(teststat)), type = "h", xlab = "chromosome number", 
            ylab = "clone statistic", pch = 18, main = paste(titles, 
                collapse = " vs "), xlim = c(0, max(cumsum(chrominfo$length), 
                kb.loc, rm.na = TRUE)), xaxt = "n")
        axis(side = 1, at = kb.loc[ind][1], label = "", tick = FALSE)
        st.now <- rev(st)
        pal.now <- rev(col.scheme$pal)
        if (length(st.now) > 0) 
            abline(h = st.now, col = pal.now, lty = 2)
        abline(v = cumsum(chrominfo$length), col = col.scheme$abline3)
        abline(v = chrom.centr, lty = 2, col = col.scheme$abline4)
        for (i in seq(1, numaut, b = 2)) mtext(paste("", i), 
            side = 1, at = chrom.mid[i], line = 0.3, col = col.scheme$mtext, 
            cex.main = 0.5)
        for (i in seq(2, numaut, b = 2)) mtext(paste("", i), 
            side = 3, at = chrom.mid[i], line = 0.3, col = col.scheme$mtext, 
            cex.main = 0.5)
        if (X) 
            if (is.even(numaut)) 
                mtext("X", side = 1, at = (chrom.mid[numaut + 
                  1]), line = 0.3, col = col.scheme$mtext, cex.main = 0.5)
            else mtext("X", side = 3, at = (chrom.mid[numaut + 
                1]), line = 0.3, col = col.scheme$mtext, cex.main = 0.5)
        if (Y) 
            if (is.even(numaut)) 
                mtext("Y", side = 3, at = (chrom.mid[numaut + 
                  2]), line = 0.3, col = col.scheme$mtext, cex.main = 0.5)
            else mtext("Y", side = 1, at = (chrom.mid[numaut + 
                2]), line = 0.3, col = col.scheme$mtext, cex.main = 0.5)
    }
}
