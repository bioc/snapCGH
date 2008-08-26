"prop.na" <- function(x)
    mean(is.na(x))



"processCGH" <- function (input, maxChromThreshold = 22, minChromThreshold = 1, 
    method.of.averaging = NULL, ID = "ID",prop.missing=0.1) 
{
    if (is.null(input$design)) {
        stop("$design component is null")
    }
    if (class(input) == "RGList") {
        MA = MA.RG(input)
        MA$design = input$design
    }
    else if (class(input) == "MAList") 
        MA = input
    else stop("Input must be of class RGList or MAList")
    for (i in 1:length(MA$design)) {
        temp <- MA$design[i] * MA$M[, i]
        MA$M[, i] <- temp
    }
    ord <- order(MA$genes$Chr, MA$genes$Position)
    if (length(which(colnames(MA$genes) == ID)) == 0) {
        stop("Specified ID column in $genes does not exist. Please check the ID argument")
    }
    colnames(MA$genes)[which(colnames(MA$genes) == ID)] = "ID"
    if (!is.null(MA$genes$Status)) {
        valStore <- attr(MA$genes$Status, "values")
        colStore <- attr(MA$genes$Status, "col")
    }
    MA <- MA[ord, ]
    ind.unmap <- which(MA$genes$Chr > maxChromThreshold | is.na(MA$genes$Chr) | 
        is.na(MA$genes$Position) | MA$genes$Chr < minChromThreshold)
    if (length(ind.unmap) > 0) {
        MA <- MA[-ind.unmap, ]
    }
    prop.miss <- apply(MA$M, 1, prop.na)
    MA <- MA[prop.miss < prop.missing, ]
    tbl <- table(MA$genes$ID)
    if (!is.null(method.of.averaging)){
      if (any(tbl > 1)) {
        tbl <- tbl[tbl > 1]
        nms <- names(tbl)
        cat("\nAveraging duplicated clones\n")
        for (i in 1:length(tbl)) {
          ind1 <- which(MA$genes$ID == nms[i])
          vec <- apply(as.matrix(MA$M[ind1, ]), 2, method.of.averaging, 
                       na.rm = TRUE)
          for (j in 1:length(ind1)) {
            if (ncol(log2ratios(MA)) > 1) {
              MA$M[ind1[j], ] <- vec
            }
            else {
              MA$M[ind1[j]] <- vec
            }
          }
        }
        dupl <- duplicated(MA$genes$ID)
        segList <- new("SegList")
        segList$M.observed <- MA$M[!dupl, , drop = FALSE]
        segList$genes <- MA$genes[!dupl, ]
      }
    }
    else {
      segList <- new("SegList")
      segList$M.observed <- MA$M
      segList$genes <- MA$genes
    }
    rownames(segList$genes) <- c(1:length(segList$genes$ID))
    if (!is.null(segList$genes$Status)) {
        attr(segList$genes$Status, "values") <- valStore
        attr(segList$genes$Status, "col") <- colStore
    }
    if (!is.null(method.of.averaging)) {
        seg.imputed <- imputeMissingValues(segList, chrominfo = chrominfo.Mb, 
            maxChrom = maxChromThreshold, smooth = 0.1)
        segList$M.observed <- seg.imputed$M
    }
    segList$design <- MA$design
    segList
}




    
"imputeMissingValues" <- 
function (seg, chrominfo = chrominfo.Mb, maxChrom = 23, 
    smooth = 0.1) 
{
    data.imp <- log2ratios <- log2ratios(seg)
    clones.info <- seg$genes
    uniq.chrom <- unique(clones.info$Chr)
    for (j in uniq.chrom[uniq.chrom <= maxChrom]) {
        cat("Processing chromosome ", j, "\n")
        centr <- chrominfo$centromere[j]
        indl <- which(clones.info$Chr == j & clones.info$Position < 
            centr)
        indr <- which(clones.info$Chr == j & clones.info$Position > 
            centr)
        kbl <- clones.info$Position[indl]
        kbr <- clones.info$Position[indr]
        for (i in 1:ncol(log2ratios)) {
            if (length(indl) > 0) {
                vecl <- log2ratios[indl, i]
                if (length(vecl[!is.na(vecl) == TRUE])!= 0)  ind <- which(!is.na(vecl)) else {ind <- 0}
                if (length(ind) > 2){ 
                  data.imp[indl, i][-ind] <- approx(lowess(kbl[ind], 
                    vecl[ind], f = smooth), xout = kbl[-ind])$y
              }}
            if (length(indr) > 0) {
               vecr <- log2ratios[indr, i]
               if (length(vecr[!is.na(vecr) == TRUE])!= 0)   ind <- which(!is.na(vecr)) else {ind <- 0}
                if (length(ind) > 2){ 
                  data.imp[indr, i][-ind] <- approx(lowess(kbr[ind], 
                    vecr[ind], f = smooth), xout = kbr[-ind])$y
            }}
        }
    }
    prop.miss <- apply(data.imp, 2, prop.na)
    if (max(prop.miss, na.rm = TRUE) > 0) {
        for (i in 1:ncol(data.imp)) {
            vec <- data.imp[, i]
            ind <- which(is.na(vec))
            if (length(ind) > 0) {
                vec[ind] <- sapply(ind, function(i) {
                  chr <- clones.info$Chr[i]
                  kb <- clones.info$Position[i]
                  if (kb >= chrominfo$centromere[chr]) 
                    median(vec[clones.info$Chr == chr & clones.info$Position >= 
                      chrominfo$centromere[chr]], na.rm = TRUE)
                  else median(vec[clones.info$Chr == chr & 
                    clones.info$Position < chrominfo$centromere[chr]], 
                    na.rm = TRUE)
                })
                vec[is.na(vec)] <- 0
                data.imp[, i] <- vec
            }
        }
    }
    prop.miss <- apply(data.imp, 2, prop.na)
    if (max(prop.miss) > 0) 
        print(paste("Missing values still remain in samples ", 
            which(prop.miss > 0)))
    seg$M <- as.matrix(data.imp)
    seg
  }


