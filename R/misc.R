##"RG.intensity.correction" <-
##function (RG, control.ID = "Dros", method.of.averaging = median, 
##    multiplier = 2, background.method = "minimum") 
##{
##  no.samples <- ncol(RG$R)
##  RG1a <- backgroundCorrect(RG, method = background.method)
##    RG1b <- RG1a
##    medR <- vector()
##    medG <- vector()
##    control.ind <- grep(control.ID, RG$genes$Name, value = TRUE)
##    control.ind1 <- unique(control.ind)
##    for (i in 1:length(no.samples)) {
##        medR[i] <- method.of.averaging(RG1a$R[RG1a$genes$Name %in% 
##            control.ind1, i])
##    }
##    for (i in 1:length(no.samples)) {
##        medG[i] <- method.of.averaging(RG1a$G[RG1a$genes$Name %in% 
##            control.ind1, i])
##   }
##    for (i in 1:length(no.samples)) {
##        RG1a$R[, i] <- as.numeric(ifelse(RG1a$R[, i] < multiplier * 
##            medR[i], ifelse(RG1a$G[, i] < multiplier * medG[i], 
##            "NA", RG1a$R[, i]), RG1a$R[, i]))
##    }
##    for (i in 1:length(no.samples)) {
##        RG1b$G[, i] <- as.numeric(ifelse(RG1b$G[, i] < multiplier * 
##            medG[i], ifelse(RG1b$R[, i] < multiplier * medR[i], 
##            "NA", RG1b$G[, i]), RG1b$G[, i]))
##    }
##    RG1a$G <- RG1b$G
##    RG1a
##}


##"MA.intensity.correction" <-
##function (current.MA = MA.hopeful, control.ID = "dros", method.of.averaging = median, 
##    multiplier = 1) 
##{
##    RG.alpha <- RG.MA(current.MA)
##    RG.beta <- RG.MA(current.MA)
##    medR <- vector()
##    medG <- vector()
##    control.ind <- grep(control.ID, RG.alpha$genes$Name, value = TRUE)
##    control.ind1 <- unique(control.ind)
##    for (i in 1:ncol(RG.alpha$R)) {
##        medR[i] <- method.of.averaging(RG.alpha$R[RG.alpha$genes$Name %in% 
##            control.ind1, i], na.rm = TRUE)
##    }
##    for (i in 1:ncol(RG.alpha$R)) {
##        medG[i] <- method.of.averaging(RG.alpha$G[RG.alpha$genes$Name %in% 
##            control.ind1, i], na.rm = TRUE)
##    }
##    for (i in 1:ncol(RG.alpha$R)) {
##        RG.alpha$R[, i] <- as.numeric(ifelse(RG.alpha$R[, i] < 
##            multiplier * medR[i], ifelse(RG.alpha$G[, i] < multiplier * 
##            medG[i], "NA", RG.alpha$R[, i]), RG.alpha$R[, i]))
##    }
##    for (i in 1:ncol(RG.alpha$R)) {
##        RG.beta$G[, i] <- as.numeric(ifelse(RG.beta$G[, i] < 
##            multiplier * medG[i], ifelse(RG.beta$R[, i] < multiplier * 
##            medR[i], "NA", RG.beta$G[, i]), RG.beta$G[, i]))
##    }
##    RG.alpha$G <- RG.beta$G
##    MA.new <- MA.RG(RG.alpha, log.transform = FALSE, bc.method = "none", 
##        offset = 0)
##    MA.new
##}

"read.clonesinfo" <-
function (file, RG, path = NULL, sep="\t", quote="\"") 
{
    if (!is.null(path)) 
        file <- file.path(path, file)
    Table <- read.table(file, sep = sep, header = T, quote=quote, as.is=TRUE, fill=TRUE)
    Chr <- as.integer(as.character(Table$Chr))
    Position <- Table$Position
    RG$genes <- data.frame(RG$genes, Position, Chr)
    RG
}

"mu1.func" <- 
function (p) 
{
    matrix(p, nrow = 1)
}

"summarizeClones" <- 
function (MA, resT = NULL, pheno = rep(1, ncol(MA)), 
    rsp.uniq = unique(pheno), thres = 0.25, factor = 2.5, all = length(rsp.uniq) == 
        1 && is.null(resT), titles = if (all) "all" else rsp.uniq) 
{

#sd.samples takes an HMMList, the summarize function takes an aCGHList.
#Can either change the summarize function to take both or ignore sd.samples

#  if (!is.null(sd.samples(aCGH.obj))) {
#        thres <- factor * (sd.samples(aCGH.obj)$madGenome)
#    }

  data <- log2ratios(MA)
    datainfo <- MA$genes
    rsp.uniq <- sort(rsp.uniq)
    colmatr <- if (length(rsp.uniq) > 1) 
        t(sapply(rsp.uniq, function(rsp.uniq.level) ifelse(pheno == 
            rsp.uniq.level, 1, 0)))
    else matrix(rep(1, length(pheno)), ncol = length(pheno), 
        nrow = 1)
    data.thres <- as.matrix(threshold(data, thresAbs = thres))
    bac.summary <- table.bac.func(dat = data.thres, colMatr = colmatr)
    if (!is.null(resT)) {
        res <- resT[order(resT$index), ]
        bac.summary <- cbind(bac.summary, res$teststat, res$rawp, 
            res$adjp)
    }
    bac.summary <- as.data.frame(bac.summary)
    nms <- c("NumPresent", "NumGain", "NumLost", "PropPresent", 
        "PropGain", "PropLost")
    cnames <- colnames(bac.summary)
    cnames[1:6] <- paste(nms, "All", sep = ".")
    if (nrow(colmatr) > 1) 
        for (m in 1:length(rsp.uniq)) cnames[(6 * m + 1):(6 * 
            (m + 1))] <- paste(nms, titles[m], sep = ".")
    if (!is.null(resT)) 
        cnames[(ncol(bac.summary) - 2):ncol(bac.summary)] <- c("stat", 
            "rawp", "adjp")
    colnames(bac.summary) <- cnames
    bac.summary <- cbind(datainfo, bac.summary)
    invisible(bac.summary)
}

"threshold" <- 
function (dat, thresAbs) 
{
    out <- matrix(0, nrow = nrow(dat), ncol = ncol(dat))
    if (length(thresAbs) == 1) 
        thresAbs <- rep(thresAbs, ncol(dat))
    if (length(thresAbs) != ncol(dat)) 
        stop("Error: number of threshold is not the same as number of\nsamples")
    sapply(1:ncol(dat), function(i) {
        tmp <- rep(0, ncol(dat))
        na.col <- is.na(dat[, i])
        col <- dat[, i][!na.col]
        tmp[!na.col] <- ifelse(col > thresAbs[i], 1, ifelse(col < 
            thresAbs[i], -1, 0))
        tmp
    })
}

log2ratios <- function(x) {
  if(!is.null(x$M.observed)){
    matrix(as.matrix(x$M.observed), nrow = nrow(x$M.observed), ncol = ncol(x$M.observed), byrow = FALSE, dimnames = dimnames(x$M.observed))
  }
  else{
    matrix(as.matrix(x$M), nrow = nrow(x$M), ncol = ncol(x$M), byrow = FALSE, dimnames = dimnames(x$M))
  }
}


