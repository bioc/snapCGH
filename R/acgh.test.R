create.resT <-
    function(resT.raw, p.adjust.method = "fdr")
{
    
    rawp <- resT.raw[ 2, ]
    adjp <- p.adjust(rawp, p.adjust.method)
    teststat <- resT.raw[ 1, ]
    
    data.frame(index = 1:ncol(resT.raw),
               teststat = teststat,
               rawp = rawp,
               adjp = adjp
               )[ order(adjp, rawp, teststat), ]

}

aCGH.test <- function(MA, rsp, test = c("survdiff", "coxph", "linear.regression"),p.adjust.method = "fdr", subset = NULL, strt = NULL, ...)
{

    l2r <- as.matrix(log2.ratios(MA))
    if (!is.null(subset))
        l2r <- l2r[ subset, ]
    test <- match.arg(test)
    pheno <- phenotype(MA)
    resT <- 
        sapply(1:nrow(l2r),
               function(i) {
                   
###                   if (i %% 100 == 0)
###                       print(i)
                   clone <- l2r[ i, ]
                   fmla <-
                       if (!is.null(strt))
                           rsp ~ clone + strata(strt)
                       else
                           rsp ~ clone
                   switch(test,
                          survdiff = {
                              
                              survdiff.fit <- try(survdiff(fmla, ...))
                              if (inherits(survdiff.fit, "try-error"))
                                  c(0, 1)
                              else
                              {
                                  
                                  etmp <- 
                                      if (is.matrix(survdiff.fit$obs))
                                          apply(survdiff.fit$exp,
                                                1,
                                                sum)
                                      else
                                          survdiff.fit$exp
                                  c(survdiff.fit$chisq,
                                    1 - pchisq(survdiff$chisq,
                                               sum(etmp > 0) - 1))
                                  
                              }
                              
                          },
                          coxph = {

                              coxph.fit <- try(coxph(fmla, ...))
                              if (inherits(coxph.fit, "try-error"))
                                  c(0, 1)
                              else
                              {
                                  cf <-  coxph.fit$coef
				  cf.se <- sqrt(coxph.fit$var)
				  cf.std <- cf/cf.se
				  c(cf.std, 2*(1-pnorm(abs(cf.std))))
                                  #logtest <-
                                  #    -2 * (coxph.fit$loglik[1] -
                                  #          coxph.fit$loglik[2])
                                  #beta <- coxph.fit$coef
                                  #df <- length(beta[!is.na(beta)])
                                  #c(logtest, 1 - pchisq(logtest, df))
                                  
                              }
                              
                          },
                          linear.regression = {
                              
				reg <- lm(fmla, ...)
				cf <- (summary(reg))$coef
				c(cf[2,3], cf[2,4])
                        ##      fstat <-
                        ##       summary(lm(clone ~ rsp, ...))$fstatistic
                        ##   c(fstat[1],
                        ##    1 - pf(fstat[1], fstat[2], fstat[3]))
                              
                          }
###                          logistic.regression = {
###
###                              glm.fit <-
###                                  try(glm(fmla, family=binomial()))
###                              if (inherits(glm.fit, "try-error"))
###                                  c(0, 1)
###                              else
###                              {
###                                  
###                                  stat <-
###                                      2 * (glm.fit$null.deviance -
###                                           glm.fit$deviance)
###                                  c(stat,
###                                    1 - pchisq(stat,
###                                               glm.fit$df.null -
###                                               glm.fit$df.residual))
###
###                              }
###                              
###                          }
                          )
                   
               }
               )
    
    create.resT(resT, p.adjust.method)
    
}

fga.func <- function(input, thres=0.25, factor=2.5, samplenames = colnames(input$M), chrominfo=chrominfo.basepair) {

#        #check if sd.samples are non-empty:
	#if (!is.null(sd.samples(aCGH.obj)) && (factor > 0)) {
#		thres <- factor*(sd.samples(aCGH.obj)$madGenome)
#	}
    
	data <- log2.ratios(input)
	data.thres <- threshold.func(data, thresAbs=thres)
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

gainLoss <-
    function (input, cols = c(1:ncol(input)), thres = 0.25)
{

    if (length(thres) == 1)
        thres <- rep(thres, ncol(input))
    if (length(thres) != ncol(input))
        stop("Error: number of thresholds is not the same as number\
of samples")
    if(class(input) == "aCGHList"){
      dt <- input$M
    }
    else if (class(input) == "SegList"){
      dt <- input$M.predicted
    }
    else{
      dt <- as.matrix(input[,cols])
    }
#    dt <- as.matrix(input[ ,cols ])
#    dt <- matrix(input[,cols], nrow = nrow(input), ncol= length(cols), byrow = FALSE)
    thr <- thres[cols]
    loss <- gain <- rep(0, nrow(dt))
    
    for (i in 1:nrow(dt))
        if (!all(is.na(dt[ i, ])))
        {
            
            x <- dt[ i, ]
            th <- thr[!is.na(x)]
            x <- x[!is.na(x)]
            tmp.gain <- x >= th
            gain[i] <- mean(tmp.gain)
            #if (any(tmp.gain))
            #    gain.med[i] <- quantile(x[tmp.gain], 1 - quant)
            tmp.loss <- x <= -th
            loss[i] <- mean(tmp.loss)
            #if (any(tmp.loss))
            #    loss.med[i] <- quantile(x[tmp.loss], quant)
            
        }
    
    list(gainP = gain,lossP = loss)    
}

"plotFreqStat" <-
function (segList, resT = NULL, pheno = rep(1, ncol(log2.ratios(segList))), 
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
    data <- log2.ratios(segList)
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

is.even <-
function (x) 
{
    if (is.numeric(x)) {
        if (x%%2 == 0) {
            TRUE
        }
        else {
            FALSE
        }
    }
    else {
        print("Warning: Input must be an integer value")
    }
}
