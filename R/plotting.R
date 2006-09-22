plotSegmentedGenome <- function (..., array = 1, naut = 22, Y = FALSE, X = FALSE, status, 
    values, pch, cex, col, chrominfo = chrominfo.Mb, ylim = c(-2, 
        2), ylb = "Log2Ratio", chrom.to.plot = NA, xlim = c(0, 
        NA), colors = NULL, mark.regions = FALSE, main = NA) 
{  
  objects <- list(...)
  nobjects <- length(objects)
  for (i in 1:nobjects) {
    if (class(objects[[i]]) != "SegList") {
      stop("Objects must be of class SegList")
    }
  }
  if (is.null(colors)) {
    colors = rep(c("blue"), nobjects)
  }
  genomePlot(input = objects[[1]], array = array, naut = naut, 
             Y = Y, X = X, status = status, values = values, chrominfo = chrominfo, 
             ylim = ylim, ylb = ylb, chrom.to.plot = chrom.to.plot, 
             xlim = xlim, pch = pch, col = col, cex = cex, main = main)
  nchr <- naut
  if (X) 
    nchr <- nchr + 1
  if (Y) 
    nchr <- nchr + 1
  chrom.start <- c(0, cumsum(chrominfo$length))[1:nchr]
  if (!is.na(chrom.to.plot)) {
    for (k in 1:nobjects) {
      current <- objects[[k]][objects[[k]]$genes$Chr == 
                              chrom.to.plot, ]
      breakpoints <- findBreakPoints(current, array)
      dup.breaks <- breakpoints[duplicated(breakpoints)]
      out <- rep(1, length(current$M.predicted[, array]) - 
                 1)
      out[dup.breaks - 1] <- 3
      out[dup.breaks] <- 3
      width <- rep(2, length(current$M.predicted[, array]) - 
                   1)
      width[dup.breaks - 1] <- 0.5
      width[dup.breaks] <- 0.5
      for (j in 2:length(objects[[k]][objects[[k]]$genes$Chr == 
                                      chrom.to.plot, array])) {
        segments(y0 = current$M.predicted[(j - 1), array], 
                 x0 = (current$genes$Position[(j - 1)]),
                 y1 = current$M.predicted[(j), array],
                 x1 = (current$genes$Position[(j)]), 
                 col = colors[k], lwd = width[j - 1], lty = out[j - 
                                                        1])
      }
      points(current$genes$Position[dup.breaks], current$M.observed[dup.breaks, 
                                                                    array], col = colors[k], pch = 16, cex = 1)
    }
  }
  else {
    for (k in 1:nobjects) {
      for (i in 1:nchr) {
        current <- objects[[k]]
        breakpoints <- findBreakPoints(current, array)
        breakpoints <- breakpoints[which(current$genes$Chr[breakpoints] == i)]
        
        segments(y0 = current$M.predicted[breakpoints[-length(breakpoints)], array],
                 x0 = (current$genes$Position[breakpoints[-length(breakpoints)]] + chrom.start[i]),
                 y1 = current$M.predicted[breakpoints[-1], array],
                 x1 = (current$genes$Position[breakpoints[-1]] + chrom.start[i]),
                 col = colors[k], lwd = 2)
      }
    }
    if(mark.regions == TRUE) {
      for(k in 1:nobjects) {
        if(!is.null(objects[[k]]$regions))
          regions <- objects[[k]]$regions
        else
          break
        for(i in 1:nchr) {
          current <- objects[[k]]
          breakpoints <- findBreakPoints(current, array)
          breakpoints <- breakpoints[which(current$genes$Chr[breakpoints] == i)]
          br.start <- breakpoints[(1:(length(breakpoints)/2)) * 2 -1]
          br.end <- breakpoints[(1:(length(breakpoints)/2)) * 2]
          
          current.regions <- regions[[array]][regions[[array]]$chr == i,]
          current.regions <- current.regions[!is.na(current.regions$color),]
          new.start <- new.end <- new.color <- NULL
          
          for(b in 1:length(br.start))
            for(r in 1:nrow(current.regions))
              {
                overlap <- intersect(current.regions$region.start[r]:current.regions$region.end[r], br.start[b]:br.end[b])
                if(length(overlap) > 0)
                  {
                    new.start <- c(new.start, min(overlap))
                    new.end <- c(new.end, max(overlap))
                    new.color <- c(new.color, current.regions$color[r])
                  }
              }
          
          if(!is.null(new.start))
            segments(y0 = current$M.predicted[new.start, array],
                     x0 = (current$genes$Position[new.start] + chrom.start[i]),
                     y1 = current$M.predicted[new.end, array],
                     x1 = (current$genes$Position[new.end] + chrom.start[i]),
                     col = new.color, lwd = 2)
        }
      }
    }   
  }
}


"genomePlot" <-
function (input, array = 1, naut = 22, 
    Y = FALSE, X = FALSE, status, values, pch, cex, col, chrominfo = chrominfo.Mb, 
    ylim = c(-2, 2), ylb = "Log2Ratio", chrom.to.plot = NA, xlim=c(0,NA), ...) 
{

  ##MALists haven't been adjust with respect to which channel is the test.
  ##This is done here, but isn't neccessary for SegLists
  if(class(input) == "MAList"){
    if (is.null(input$design)) 
        stop("MA$design component is null")
    for(i in 1:length(input$design)){
      temp <- input$design[i]* input$M[,i]
      input$M[,i] <- temp
    }
  }
  else if(class(input) == "SegList"){} 
  else{
    stop("Class must be either MAList or SegList")
  }
   
  
  data <- log2ratios(input)
  datainfo <- input$genes
    
  ord <- order(datainfo$Chr, datainfo$Position)
  chrom <- datainfo$Chr[ord]
  kb <- datainfo$Position[ord]
  name <- (colnames(data))[array]

  data <- matrix(data[ord, ], nrow = nrow(as.matrix(data[ord,])),ncol = ncol(data), b = FALSE, dimnames = dimnames(data))
  ind.unmap <- which(chrom < 1 | is.na(chrom) | is.na(kb) | (as.numeric(chrom) > (naut + 2)))
  
  if (missing(status)) status <- input$genes$Status
  
  if (length(ind.unmap) > 0) {
        chrom <- chrom[-ind.unmap]
        kb <- kb[-ind.unmap]
        data <- as.matrix(data[-ind.unmap,])
        ## code dealing with the spot types functionality
        valStore <- attr(status,"values")
        colStore <- attr(status,"col")
        status <- status[ord]
        status <- status[-ind.unmap]
        attr(status,"values") <- valStore
        attr(status,"col") <- colStore
      }

  nchr <- naut
  if (X) nchr <- nchr + 1
  if (Y) nchr <- nchr + 1
  
  if (!is.na(chrom.to.plot)){
    status <-  NULL
    nchr <- chrom.to.plot 
    data <- matrix(data[chrom==nchr,], nrow = nrow(as.matrix(data[chrom==nchr,])), ncol = ncol(data), b = FALSE, dimnames = dimnames(data[chrom==nchr,]))
    clone.genomepos <- kb[chrom == nchr]
    chrom <- chrom[chrom == nchr]
    chrominfo <- chrominfo[nchr, ]
    chrominfo <- chrominfo[1:nchr, ]
    chrom.start <- 0
    par(xaxt = "s", cex = 0.6, pch = 18, lab = c(6, 6, 7), cex.axis = 1.5, xaxs = "i")
    } else {
      data <- matrix(data[chrom <= nchr, ], nrow = nrow(as.matrix(data[chrom <= nchr,])),ncol = ncol(data), byrow = FALSE, dimnames = dimnames(data[chrom <= nchr,])) 
      kb <- kb[chrom <= nchr]
      chrom <- chrom[chrom <= nchr]
      chrominfo <- chrominfo[1:nchr, ]
      chrom.start <- c(0, cumsum(chrominfo$length))[1:nchr]
      clone.genomepos <- vector()
      for (i in 1:nchr) {clone.genomepos[chrom == i] <- kb[chrom == i] + chrom.start[i]}
      par(xaxt = "n", tck = -0, cex = 0.6, pch = 18, lab = c(1, 6, 7), cex.axis = 1.5, xaxs = "i")
    }

  chrom.centr <- chrom.start + chrominfo$centr
  chrom.mid <- chrom.start + chrominfo$length/2
  x <- clone.genomepos
  y <- data[, array]
  if (is.na(xlim[2])) {xlim[2] <- clone.genomepos[sum(clone.genomepos > 0)]}
  
   #Plotting functions 
 
  if (!is.na(chrom.to.plot)){
    plot(x, y, ylim = ylim, xlim = xlim, xlab = "Distance along chromosome (Mb)", ylab = "" , col = "black", bg="white")
    mtext(chrom.to.plot, side = 1, line = 0.3, col = "red")}  else {
      plot(x, y, ylim = ylim, xlab = "", ylab = "", xlim = xlim, col = "black", bg="white")
      for (i in seq(1, naut, b = 2)) mtext(i, side = 1, at = chrom.mid[i], line = 0.3, col = "red")
      for (i in seq(2, naut, b = 2)) mtext(i, side = 3, at = chrom.mid[i], line = 0.3, col = "red")}

  title(main = paste(array, " ", name), ylab = ylb, xlab = "", cex.lab = 1.5)    
  if (X & is.na(chrom.to.plot)){ mtext("X", side = 1, at = chrom.mid[naut + 1], line = 0.3, col = "red")}
  if (Y & is.na(chrom.to.plot)){ mtext("Y", side = 3, at = chrom.mid[naut + 2], line = 0.3, col = "red")}


  
  abline(v = c(chrom.start, (chrom.start[nchr] + chrominfo$length[nchr])), lty = 1)
  abline(h = seq(min(ylim), max(ylim), b = 0.5), lty = 3)
  abline(v = (chrominfo$centromere + chrom.start), lty = 3, col = "red")

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

#I can't get the dendrogram section of this to work.
#The matrix transpose screws it completely as the dist function
#returns a single value and the plotting function doesn't accept the
#hcl object.

"heatmapGenome" <-
function (input, response = as.factor(rep("All", ncol(input))), 
    chrominfo = chrominfo.Mb, cutoff = 1, lowCol = "blue", 
   highCol = "yellow", midCol = "white", ncolors = 50, byclass = FALSE, 
    showaber = FALSE, amplif = 1, homdel = -0.75, samplenames = colnames(input), 
    vecchrom = 1:22, titles = "Image Plot", methodS = "ward", 
    categoricalPheno = TRUE, CENTROMERE = FALSE) 
{

  ##MALists haven't been adjust with respect to which channel is the test.
  ##This is done here, but isn't neccessary for SegLists
  if(class(input) == "MAList"){
    if (is.null(input$design)) 
        stop("MA$design component is null")
    for(i in 1:length(input$design)){
      temp <- input$design[i]* input$M[,i]
      input$M[,i] <- temp
    }
	data <- input$M
  }
  else if(class(input) == "SegList"){
	data <- input$M.predicted
	} 
  else{
    stop("Class must be either MAList or SegList")
  }
    
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
 #       data <- log2ratios.imputed(input)
 #   else data <- log2ratios(input)
 #   data <- log2ratios(input)
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
##    genomePlot(input, samples = samples, naut = naut, 
##        Y = Y, X = X, status = status, values = values, pch = pch, col = col, cex = cex, chrominfo = chrominfo, chrom.to.plot = chrom.to.plot, 
##        Z = Z, xlower = xlower, xupper = xupper)
##    identify(input$genes[input$genes$Chr == chrom.to.plot, colnames(input$genes) == "Position"]/1000,
##             log2ratios(input)[input$genes$Chr == chrom.to.plot, samples], labels = label)
##}

