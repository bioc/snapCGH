
"processCGH" <-
function (MA, chrom.remove.threshold = 22, chrom.below.threshold = 1, method.of.averaging = NULL, ID = "ID") 
{
    if (is.null(MA$design)) 
        stop("MA object doesn't contain design component")

    ord <- order(MA$genes$Chr, MA$genes$Position) # re-ordering the clones by chromosome and position on a chromosome
    colnames(MA$genes)[which(colnames(MA$genes) == ID)] = "ID" #renaming the the specified column to "ID"

#Code to stop the Status attributes being lost when the MAList is reordered.
    if(!is.null(MA$genes$Status)){
      valStore <- attr(MA$genes$Status, "values")
      colStore <- attr(MA$genes$Status, "col")
    }
      MA <- MA[ord,]
    
    ind.unmap <- which(MA$genes$Chr > chrom.remove.threshold | is.na(MA$genes$Chr) | is.na(MA$genes$Position) | MA$genes$Chr < chrom.below.threshold)
        if (length(ind.unmap) > 0) {
          MA <- MA[-ind.unmap, ]
          }
    prop.miss <- apply(MA$M, 1, prop.na)
    MA <- MA[prop.miss < 0.1,]
    ## Removing duplicated clones

    tbl <- table(MA$genes$ID)
    if (any(tbl > 1)) {
        tbl <- tbl[tbl > 1]
        nms <- names(tbl)
        if (!is.null(method.of.averaging)) {
          cat("\nAveraging duplicated clones\n")
          for (i in 1:length(tbl)) {
            ind1 <- which(MA$genes$ID == nms[i])
            vec <- apply(as.matrix(MA$M[ind1,]), 2, method.of.averaging, na.rm = TRUE)
            if(!is.null(MA$weights)){
              vec2 <-  apply(as.matrix(MA$M[ind1,]), 2, method.of.averaging, na.rm = TRUE)
            }
            for (j in 1:length(ind1)) {
              if (ncol(log2.ratios(MA)) > 1) {
                MA$M[ind1[j], ] <- vec
                if(!is.null(MA$weights)){
                  MA$weights[ind1[j], ] <- vec2
                }
                else {
                  MA$M[ind1[j]] <- vec
                  if(!is.null(MA$weights)){
                    MA$weights[ind1[j]] <- vec2
                  }
                }
              }
            }
            dupl <- duplicated(MA$genes$ID)
            MA$genes <- MA$genes[!dupl, ]
            MA$M <- MA$M[!dupl, ,drop = FALSE]
            if(!is.null(MA$weights)){
              MA$weights <- MA$weights[!dupl,,drop = FALSE]
            }
          }
        }
        MA$printer <- NULL
        MA$A <- NULL

	MA$genes$ID <- factor(MA$genes$ID)
    	rownames(MA$genes) <- c(1:length(MA$genes$ID))
    	if(!is.null(MA$genes$Status)){
            attr(MA$genes$Status, "values") <- valStore
            attr(MA$genes$Status, "col") <- colStore
          }
      }
                                        # The imputation step

    if (!is.null(method.of.averaging)){
    MA.imputed <- impute.lowess(MA, chrominfo = chrominfo.basepair, maxChrom = chrom.remove.threshold, smooth = 0.1)
    MA$M <- MA.imputed$M}
    MA

  }
  

"impute.lowess" <- 
function (MA, chrominfo = chrominfo.basepair, maxChrom = 23, 
    smooth = 0.1) 
{
    data.imp <- log2.ratios <- log2.ratios(MA)
    clones.info <- MA$genes
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
        for (i in 1:ncol(log2.ratios)) {
            if (length(indl) > 0) {
                vecl <- log2.ratios[indl, i]
                if (length(vecl[!is.na(vecl) == TRUE])!= 0)  ind <- which(!is.na(vecl)) else {ind <- 0}
                if (length(ind) > 2){ 
                  data.imp[indl, i][-ind] <- approx(lowess(kbl[ind], 
                    vecl[ind], f = smooth), xout = kbl[-ind])$y
              }}
            if (length(indr) > 0) {
               vecr <- log2.ratios[indr, i]
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
    MA$M <- as.matrix(data.imp)
    MA
  }

#selection of little accessor methods.
#Might put these somewhere else in a bit

lengthGain.na <-
    function(x)
    sum(x == 1, na.rm = TRUE)

propGain.na <-
    function(x)
    mean(x == 1, na.rm = TRUE)

lengthLoss.na <-
    function(x)
    sum(x == -1, na.rm = TRUE)

propLoss.na <-
    function(x)
    mean(x == -1, na.rm = TRUE)

prop.na <-
    function(x)
    mean(is.na(x))

