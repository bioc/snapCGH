"plotGenome" <-
function (input, array = 1, naut = 22, 
    Y = TRUE, X = TRUE, status, values, pch, col, cex, chrominfo = chrominfo.basepair, 
    yScale = c(-2, 2), samplenames = (colnames(input))[array], 
    ylb = "Log2Ratio", title.font = 2, chrom.to.plot = 1, Z = FALSE, 
    xlower = 0, xupper = clone.genomepos[sum(clone.genomepos > 0)]/1000) 
{
  
  data = log2.ratios(input)
  datainfo <- input$genes
    nchr <- naut
    if (X) 
        nchr <- nchr + 1
    if (Y) 
        nchr <- nchr + 1
    if (Z) 
        nchr <- chrom.to.plot
    
    ord <- order(datainfo$Chr, datainfo$Position)
    chrom <- datainfo$Chr[ord]
    kb <- datainfo$Position[ord]

    data <- matrix(data[ord, ], nrow = nrow(as.matrix(data[ord,])),
                   ncol = ncol(data), byrow = FALSE, dimnames = dimnames(data))
    ind.unmap <- which(chrom < 1 | is.na(chrom) | is.na(kb) | (chrom > (naut + 
        2)))
  
  if (missing(status)) 
          status <- input$genes$Status
  
    if (length(ind.unmap) > 0) {
        chrom <- chrom[-ind.unmap]
        kb <- kb[-ind.unmap]
        data <- matrix(data[-ind.unmap, ], nrow = nrow(as.matrix(data[-ind.unmap,])),
                   ncol = ncol(data), byrow = FALSE,  dimnames = dimnames(data))
#code dealing with the spot types functionality
          valStore <- attr(status,"values")
          colStore <- attr(status,"col")
          status <- status[-ind.unmap]
          attr(status,"values") <- valStore
          attr(status,"col") <- colStore
        }

      
    if (Z) 
        data <- matrix(data[chrom == nchr, ], nrow = nrow(as.matrix(data[chrom == nchr,])),
                   ncol = ncol(data), byrow = FALSE, dimnames = dimnames(data))
    else data <- matrix(data[chrom <= nchr, ], nrow = nrow(as.matrix(data[chrom <= nchr,])),
                   ncol = ncol(data), byrow = FALSE, dimnames = dimnames(data))
   if (Z) 
        kb <- kb[chrom == nchr]
    else kb <- kb[chrom <= nchr]
    if (Z) 
        chrom <- chrom[chrom == nchr]
    else chrom <- chrom[chrom <= nchr]
    if (Z) 
        chrominfo <- chrominfo[nchr, ]
    else chrominfo <- chrominfo[1:nchr, ]
    if (Z) 
        chrom.start <- 0
    else chrom.start <- c(0, cumsum(chrominfo$length))[1:nchr]
    chrom.centr <- chrom.start + chrominfo$centr
    chrom.mid <- chrom.start + chrominfo$length/2
    chrom.rat <- chrominfo$length/max(chrominfo$length)
    if (Z) 
        par(cex = 0.6, pch = 18, lab = c(6, 6, 7), cex.axis = 1.5, 
            xaxs = "i")
    else par(cex = 0.6, pch = 18, lab = c(1, 6, 7), cex.axis = 1.5, 
        xaxs = "i")

  k <- array
  vec <- data[, array]
  name <- samplenames
  
  clone.genomepos <- rep(0, length(kb))
  if (Z) 
    clone.genomepos <- kb + chrom.start
  else
    for (i in 1:nrow(chrominfo)) clone.genomepos[chrom == i] <- kb[chrom == i] + chrom.start[i]
  y.min <- rep(yScale[1], nrow(chrominfo))
  y.max <- rep(yScale[2], nrow(chrominfo))
  for (i in 1:nrow(chrominfo)) {
    if (minna(vec[(chrom == i)]) < y.min[i]) 
      y.min[i] <- minna(vec[(chrom == i)])
    if (maxna(vec[(chrom == i)]) > y.max[i]) 
      y.max[i] <- maxna(vec[(chrom == i)])
  }
  ygenome.min <- minna(y.min)
  ygenome.max <- maxna(y.max)
  
                                        #Testing new plotting functions
  x = clone.genomepos/1000
  y = vec
  
  if (Z) 
    plot(clone.genomepos/1000, vec, ylim = yScale, xlab = "Distance along chromosome (Mb)", 
         ylab = "", xlim = c(xlower, xupper), col = "black")
  else plot(clone.genomepos/1000, vec, ylim = yScale, xlab = "", 
            ylab = "", xlim = c(0, clone.genomepos[sum(clone.genomepos > 0)]/1000), col = "black")
  title(main = paste(k, " ", name), ylab = ylb, xlab = "", cex.lab = 1.5, cex.main = title.font)
  if (Z) 
    mtext(paste("", chrom.to.plot), side = 1, at = (chrom.mid/1000), 
          line = 0.3, col = "red")
  else
    for (i in seq(1, naut, b = 2)) mtext(paste("", i), side = 1, at = chrom.mid[i]/1000, line = 0.3, col = "red")
  if (Z) ""
  else
    for (i in seq(2, naut, b = 2)) mtext(paste("", i), side = 3, at = chrom.mid[i]/1000, line = 0.3, col = "red")
  if (X) 
    mtext("X", side = 1, at = chrom.mid[naut + 1]/1000, line = 0.3, col = "red")
  if (Y) 
    mtext("Y", side = 3, at = chrom.mid[naut + 2]/1000, line = 0.3, col = "red")
  abline(v = c(chrom.start/1000, (chrom.start[nrow(chrominfo)] + 
           chrominfo$length[nrow(chrominfo)])/1000), lty = 1)
  abline(h = seq(min(yScale), max(yScale), b = 0.5), lty = 3)
  abline(v = (chrominfo$centromere + chrom.start)/1000, 
         lty = 3, col = "red")

                                        #Code stolen from limma to use the spottype functionality
  if(is.null(status) || all(is.na(status))) {
    if(missing(pch)) pch=16
    if(missing(cex)) cex=0.3
    points(x,y,pch=pch[[1]],cex=cex[1])
  } else {
    if(missing(values)) {
      if(is.null(attr(status,"values")))
                                        #				values <- names(sort(table(status),decreasing=TRUE))
        values <- as.character(status[,1])
      else
        values <- attr(status,"values")
    }
                                        #		Non-highlighted points
    sel <- !(status %in% values)
    nonhi <- any(sel)
    if(nonhi) points(x[sel],y[sel],pch=16,cex=0.3)
    
    nvalues <- length(values)
    
    if(missing(pch)) {
      if(is.null(attr(status,"pch")))
        pch <- rep(16,nvalues)
      else
        pch <- attr(status,"pch")
    }
    if(missing(cex)) {
      if(is.null(attr(status,"cex"))) {
        cex <- rep(1,nvalues)
        if(!nonhi) cex[1] <- 0.3
      } else
      cex <- attr(status,"cex")
    }
    if(missing(col)) {
      if(is.null(attr(status,"col"))) {
        col <- nonhi + 1:nvalues
      } else
      col <- attr(status,"col")
    }
    pch <- rep(pch,length=nvalues)
    col <- rep(col,length=nvalues)
    cex <- rep(cex,length=nvalues)
    for (i in 1:nvalues) {
      sel <- status==values[i]
      points(x[sel],y[sel],pch=pch[[i]],cex=cex[i],col=col[i])
    }
  } 
}


plotSummaryProfile <-
    function(input, geList,
             response = as.factor(rep("All", ncol(input))),
             titles = unique(response[!is.na(response)]), X = TRUE,
             Y = FALSE, maxChrom = 23,
             chrominfo = chrominfo.basepair,
             num.plots.per.page = length(titles), factor = 2.5, thresAbs=100)
{

    ind.samp <- which(!is.na(response))
    resp.na <- response[ind.samp]
    response.uniq <- sort(unique(resp.na))

    df.not.na <-
        data.frame(response = response,
                   numtrans =
                   apply(geList$num.transitions, 2, sum, na.rm = TRUE),
                   numtrans.binary =
                   apply(geList$num.transitions.binary, 2, sum, na.rm = TRUE),
                   numaber =
                   apply(geList$num.aberrations, 2, sum, na.rm = TRUE),
                   numaber.binary =
                   apply(geList$num.aberrations.binary, 2, sum, na.rm = TRUE),
                   numamplif =
                   apply(as.matrix(geList$num.amplifications[ 1:maxChrom, ]), 2, sum,
                         na.rm = TRUE),
                   numamplif.binary =
                   apply(as.matrix(geList$num.amplifications.binary[ 1:maxChrom, ]),
                         2, sum, na.rm = TRUE),
                   numoutlier =
                   apply(geList$num.outliers, 2, sum, na.rm = TRUE),
                   num.outliers.binary =
                   apply(geList$num.outliers.binary, 2, sum, na.rm = TRUE),
                   numchromgain =
                   apply(as.matrix(geList$whole.chrom.gain.loss[ 1:maxChrom, ]), 2,
                         length.num.func, 1),
###                   apply(ge$whole.chrom.gain[ 1:maxChrom, ], 2, sum,
###                         na.rm = TRUE),
                   numchromloss =
                   apply(as.matrix(geList$whole.chrom.gain.loss[ 1:maxChrom, ]), 2,
                         length.num.func, -1),
###                   apply(ge$whole.chrom.loss[ 1:maxChrom, ], 2, sum,
###                         na.rm = TRUE)
                   sizeamplicon =
                   apply(as.matrix(geList$size.amplicons[ 1:maxChrom, ]), 2, sum,
                         na.rm = TRUE),
                   numamplicon =
                   apply(as.matrix(geList$num.amplicons[ 1:maxChrom, ]), 2, sum,
                         na.rm = TRUE)
                   )[ which(!is.na(response)), ]
    attach(df.not.na)
    numchromchange <- numchromgain + numchromloss

    deletions <- threshold.func(log2.ratios(input), thresAbs = thresAbs)
    
    boxplot.this <-
        function(events, title, sig = 6)
        {

            p.value <-
                if (length(response.uniq) > 1)	
                    signif(kruskal.test(events ~ resp.na)$p.value,
                           sig)
                else
                    ""
            boxplot(events ~ resp.na, notch = TRUE, names = titles,
                    varwidth = TRUE, main = paste(title, p.value), cex.main=0.8)
            
        }

#############################################
    ##Plot1:
    par(mfrow = c(2, 2))

    boxplot.this(numtrans, "Number of Transitions")
    boxplot.this(numtrans.binary,
                 "Number of Chrom containing Transitions")
    boxplot.this(numaber, "Number of Aberrations")
    boxplot.this(numchromchange, "Number of Whole Chrom Changes")

#############################################
    ##Plot2:
    
    boxplot.this(numamplif, "Number of Amplifications")
    boxplot.this(numamplif.binary,
                 "Number of Chrom containing Amplifications")
    boxplot.this(numamplicon, "Number of Amplicons")
    boxplot.this(sizeamplicon, "Amount of Genome Amplified")

#############################################
    ##Plot3:

    out <- as.data.frame(fga.func(input, factor=factor,
        chrominfo=chrominfo))[ which(!is.na(response)), ]
    boxplot.this(out$gainP, "Fraction of Genome Gained")
    boxplot.this(out$lossP, "Fraction of Genome Lost")
    boxplot.this(out$gainP+out$lossP, "Fraction of Genome Altered")

#############################################

    plot.freq.this <-
        function(matr, i, ylb)
        {
            
            par(mfrow = c(num.plots.per.page, 1))

            out <-
                sapply(1:length(response.uniq),
                       function(j)
                       apply(as.matrix(matr[,which(resp.na == response.uniq[j])]),
                             1,
                             prop.num.func,
                             i
                             )
                       )
            mx <- max(c(out), na.rm = TRUE)
            if (length(titles == 1))
                out <- cbind(out, out)
            
            plotGenome(input, samples = 1:length(titles),
                       yScale = c(0, mx), naut = 22,
                       X = X, Y = Y, ylb = ylb,
                       chrominfo = chrominfo[ 1:maxChrom, ],
                       samplenames = titles
                       )
            
        }
    
    ##Plot3: trans start
    plot.freq.this(geList$transitions$trans.matrix, 1,
                  "Proportion of Transition Starts")

    ##Plot4: trans end
    plot.freq.this(geList$transitions$trans.matrix, 2,
                  "Proportion of Transition Ends")

    ##Plot5: amplification
    plot.freq.this(geList$amplifications$amplif, 1,
                  "Proportion of Amplifications")
    mtext("Amplifications", side = 3, outer = TRUE)

    ##Plot6: aberration
    plot.freq.this(geList$aberrations$aber, 1,
                   "Proportion of Aberrations")
    mtext("Aberrations", side = 3, outer = TRUE)

    ##Plot8: homozygous deletions
    plot.freq.this(deletions, -1,
                  "Proportion of Homozygous Deletions")
    mtext("Homozygous Deletions", side = 3, outer = TRUE)

    ##Plot7: whole chromosomal gain/loss:

    par(mfrow = c(num.plots.per.page, 2), lab = c(5,6,7))
    
    matr <- as.matrix(geList$whole.chrom.gain.loss[ 1:22, ])
    out.gain <- matrix(NA, nrow = nrow(matr), ncol = length(titles))
    out.loss <- matrix(NA, nrow = nrow(matr), ncol = length(titles))
    for (j in 1:length(response.uniq))
    {
        
        ind <- which(response == response.uniq[j])
        out.gain[ ,j ] <-
            apply(as.matrix(matr[ ,ind ]), 1, length.num.func, 1) / ncol(matr)
        out.loss[ ,j ] <-
            apply(as.matrix(matr[ ,ind ]), 1, length.num.func, -1) / ncol(matr)
        
    }
    mx.gain <- max(c(out.gain), na.rm = TRUE)
    mx.loss <- max(c(out.loss), na.rm = TRUE)
    mx <- max(mx.gain, mx.loss)
    for (j in 1:length(titles))
    {
        
        plot(1:22, out.gain[,j], pch = 20,
             main = as.character(titles[j]),
             xlab = "chromosome",
             ylab = "Proportion of whole chromosomes gains",
             ylim = c(0, mx), xlim = c(0,23))
        plot(1:22, out.loss[,j], pch = 20,
             main = as.character(titles[j]),
             xlab = "chromosome",
             ylab = "Proportion of whole chromosomes losses",
             ylim = c(0, mx), xlim = c(0,23))
        
    }

    detach(df.not.na)
}

#I can't get the dendrogram section of this to work.
#The matrix transpose screws it completely as the dist function
#returns a single value and the plotting function doesn't accept the
#hcl object.

"clusterGenome" <-
function (input, response = as.factor(rep("All", ncol(input))), 
    chrominfo = chrominfo.basepair, cutoff = 1, lowCol = "blue", 
   highCol = "yellow", midCol = "white", ncolors = 50, byclass = FALSE, 
    showaber = FALSE, amplif = 1, homdel = -0.75, samplenames = colnames(input), 
    vecchrom = 1:22, titles = "Image Plot", methodS = "ward", 
    #imp = TRUE, 
    categoricalPheno = TRUE, CENTROMERE = FALSE) 
{
    if(ncol(input) == 1) {
      stop("You need at least 2 samples to use this function")
    }
    input <- input[, !is.na(response)]
    samplenames <- samplenames[!is.na(response)]
    response <- response[!is.na(response)]
    if (categoricalPheno) {
        resp0 <- response
        resp0.num <- as.numeric(as.factor(resp0))
        resp <- as.numeric(as.factor(resp0))
        if (!(byclass)) {
            resp <- rep(1, length(resp0))
        }
        tbl.resp <- table(resp)
        label.col <- rainbow(length(unique(resp)))
    }
    else {
        byclass <- FALSE
        resp0 <- response
        resp0.num <- resp0
        resp <- rep(1, length(resp0))
        tbl.resp <- table(resp)
        label.col <- rainbow(length(unique(resp)))
    }
    datainfo <- input$genes
 #   if (imp) 
 #       data <- log2.ratios.imputed(input)
 #   else data <- log2.ratios(input)
    data <- log2.ratios(input)
    indUse <- vector()
    chromb <- 0
    for (i in 1:length(vecchrom)) indUse <- c(indUse, which(datainfo$Chr == 
        vecchrom[i]))
    if (CENTROMERE) 
        indUse <- indUse[datainfo$Position[datainfo$Chr == vecchrom[i]] <= 
            chrominfo$centromere[vecchrom[i]]]
    else indUse <- indUse
    for (i in 1:length(vecchrom)) if (CENTROMERE) 
        chromb <- c(chromb, length(indUse))
    else chromb <- c(chromb, length(which(datainfo$Chr == vecchrom[i])))
    chromb <- cumsum(chromb)
    datainfo <- datainfo[indUse, ]
    data <- as.matrix(data[indUse, ])
    kb <- datainfo$Position

#    if (dendPlot) {
#      cr <- dist(t(data))
#      hcl <- hclust(cr, method = methodS)
#      plot(hcl)
#    }
    dt.cp <- data
    dt <- apply(data, 2, floor.func, cutoff)
    dt <- dt[, order(resp)]
    dt.cp <- dt.cp[, order(resp)]
    resp0 <- resp0[order(resp)]
    resp0.num <- resp0.num[order(resp)]
    samplenames <- samplenames[order(resp)]
    resp <- resp[order(resp)]
    start <- 1
    ord <- rep(0, length(resp))
    for (i in 1:(length(tbl.resp))) {
        ind <- which(resp == i)
        cr <- dist(t(data[, ind]))
        ord[start:sum(tbl.resp[1:i])] <- hclust(cr, method = methodS)$ord + 
            start - 1
        start <- sum(tbl.resp[1:i]) + 1
    }
    dt <- dt[, ord]
    dt.cp <- dt.cp[, ord]
    resp <- resp[ord]
    resp0 <- resp0[ord]
    resp0.num <- resp0.num[ord]
    samplenames <- samplenames[ord]
    image(x = (1:length(kb)), y = 1:length(resp), z = dt, col = maPalette(low = lowCol, 
        high = highCol, mid = midCol, k = ncolors), axes = FALSE, 
        xlab = "", ylab = "", zlim = c(-cutoff, cutoff))
    if (showaber) {
        for (j in 1:ncol(dt)) {
            tmp <- dt.cp[, j]
            i <- (1:length(tmp))[tmp >= amplif & !is.na(tmp)]
            if (length(i) > 0) {
                points(i, rep(j, length(i)), col = "yellow", 
                  pch = 20, cex = 0.7)
            }
            i <- (1:length(tmp))[tmp <= homdel & !is.na(tmp)]
            if (length(i) > 0) {
                points(i, rep(j, length(i)), col = "skyblue", 
                  pch = 20, cex = 0.7)
            }
        }
    }
    for (j in seq(1, ncol(dt), b = 2)) {
        mtext(paste((samplenames)[j], ""), side = 2, at = j, 
            line = 0.3, col = label.col[((resp0.num)[j])], cex = 0.5, 
            las = 0)
    }
    for (j in seq(2, ncol(dt), b = 2)) {
        mtext(paste((samplenames)[j], ""), side = 4, at = j, 
            line = 0.3, col = label.col[((resp0.num)[j])], cex = 0.5, 
            las = 0)
    }
    title(xlab = "clone", ylab = "sample", main = titles)
    abline(v = chromb, col = "white", lty = 1, lwd = 1)
    loc <- chromb[-1] - diff(chromb)/2
    if (length(vecchrom) > 1) {
        for (i in seq(2, length(vecchrom), b = 2)) {
            mtext(paste("", vecchrom[i]), side = 3, at = (loc[i]), 
                line = 0.3, cex.main = 0.25)
        }
    }
    for (i in seq(1, length(vecchrom), b = 2)) {
        mtext(paste("", vecchrom[i]), side = 1, at = (loc[i]), 
            line = 0.3, cex.main = 0.25)
    }
#    if (dendPlot) {
#      print("here!")
#        if (length(unique(resp0)) > 1) {
#          plot(hcl, labels = response, main = "Dendogram")
         # plot(hcl)
#        }
#        else {
#         plot(hcl, labels = (colnames(input$M)), main = "Dendogram")
          #   plot(hcl)
#          }
#    }
}



##"labelled.genome.plot" <-
##function (input, samples = 1, naut = 22, Y = FALSE, X = FALSE,  status = NA, values, pch, col, cex, 
##    chrominfo = chrominfo.basepair, Z = TRUE, chrom.to.plot = 8, xlower = 0, 
##    xupper = max((input$genes$Position[input$genes$Chr == chrom.to.plot])/1000), identifier.to.plot = "Position") 
##{
##    label <- vector()
##    label <- input$genes[input$genes[, 
##        colnames(input$genes) == "Chr"] == chrom.to.plot, 
##       colnames(input$genes) == identifier.to.plot]
##    plotGenome(input, samples = samples, naut = naut, 
##        Y = Y, X = X, status = status, values = values, pch = pch, col = col, cex = cex, chrominfo = chrominfo, chrom.to.plot = chrom.to.plot, 
##        Z = Z, xlower = xlower, xupper = xupper)
##    identify(input$genes[input$genes$Chr == chrom.to.plot, colnames(input$genes) == "Position"]/1000,
##             log2.ratios(input)[input$genes$Chr == chrom.to.plot, samples], labels = label)
##}

"plotSegmentationStates" <-
function (segList, geList, array=1, chr = 1:length(unique(segList$genes$Chr)), 
    maxChrom = 22, chrominfo = chrominfo.basepair, 
    yScale = c(-2, 2), samplenames = sample.names(segList), 
    xlower = 0, xupper = chrominfo$length[chr[j]]/1000) 
{
    if (length(array) > 1) 
        stop("plotHmmStates currently prints only 1 sample at a\ntime\n")

    aber <- geList$aberrations$aber
    amplif <- geList$amplifications$amplif
    trans <- geList$transitions$trans.matr
    outliers <- geList$outliers$outlier
    pred <- geList$outliers$pred.out
    chrom.rat <- chrominfo$length/max(chrominfo$length)
    chrom.start <- c(0, cumsum(chrominfo$length))[1:maxChrom]
    chrom.mid <- chrom.start + chrominfo$length[1:maxChrom]/2
    chrom <- segList$genes$Chr
    par(lab = c(15, 6, 7), pch = 18, cex = 1, lwd = 1, mfrow = c(2, 
        1))
    for (j in 1:length(chr)) {
        ind.nonna <- which(!is.na(segList$M.observed[chrom == chr[j], 
            array]))
        kb <- segList$genes$Position[chrom == chr[j]][ind.nonna]/1000
        obs <- segList$M.observed[chrom == chr[j], array][ind.nonna]
        states <- segList$state[chrom == chr[j], array][ind.nonna]
        nstates <- length(unique(states))
        abernow <- aber[chrom == chr[j], array][ind.nonna]
        outliersnow <- outliers[chrom == chr[j], array][ind.nonna]
        amplifnow <- amplif[chrom == chr[j], array][ind.nonna]
        transnow <- trans[chrom == chr[j], array][ind.nonna]
        prednow <- obs
        predicted <- pred[chrom == chr[j], array][ind.nonna]
        prednow[outliersnow == 0 & abernow == 0] <- predicted[outliersnow == 
            0 & abernow == 0]
        y.min <- min(yScale[1], min(obs))
        y.max <- max(yScale[2], max(obs))
        plot(kb, obs, xlab = "", ylab = "", ylim = c(y.min, y.max), 
            type = "l", col = "blue", xlim = c(xlower, xupper))
        points(kb, obs, col = "black")
        title(main = paste("Sample", array, samplenames[array], 
            "- Chr", chr[j], "Number of states", nstates), xlab = "kb (in 1000's)", 
            ylab = "data (observed)")
        abline(h = seq(-2, 2, b = 0.5), lty = 3)
        abline(v = chrominfo$centromere[chr[j]]/1000, lty = 2, 
            col = "red", lwd = 3)
        if (nstates > 1) {
            abline(v = kb[transnow == 1], col = "blue", lwd = 2)
            abline(v = kb[transnow == 2], col = "green", lty = 2, 
                lwd = 0.5)
        }
        if (length(outliersnow[outliersnow == 1]) > 0) 
            points(kb[outliersnow == 1], obs[outliersnow == 1], 
                col = "yellow")
        if (length(abernow[abernow == 1]) > 0) 
            points(kb[abernow == 1], obs[abernow == 1], col = "orange")
        if (length(amplifnow[amplifnow == 1]) > 0) 
            points(kb[amplifnow == 1], obs[amplifnow == 1], col = "red")
        plot(kb, prednow, xlab = "", ylab = "", ylim = c(y.min, 
            y.max), type = "l", col = "blue", xlim = c(xlower, 
            xupper))
        points(kb, prednow, col = "black")
        title(xlab = "kb (in 1000's)", ylab = "data (smoothed)")
        abline(h = seq(-2, 2, b = 0.5), lty = 3)
        abline(v = chrominfo$centromere[chr[j]]/1000, lty = 2, 
            col = "red", lwd = 3)-
        if (nstates > 1) {
            abline(v = kb[transnow == 1], col = "blue", lwd = 2)
            abline(v = kb[transnow == 2], col = "green", lty = 2, 
                lwd = 0.5)
        }
        if (length(outliersnow[outliersnow == 1]) > 0) 
            points(kb[outliersnow == 1], obs[outliersnow == 1], 
                col = "yellow")
        if (length(abernow[abernow == 1]) > 0) 
            points(kb[abernow == 1], obs[abernow == 1], col = "orange")
        if (length(amplifnow[amplifnow == 1]) > 0) 
            points(kb[amplifnow == 1], obs[amplifnow == 1], col = "red")

      }
    
}

