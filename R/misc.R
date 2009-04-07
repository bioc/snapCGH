read.clonesinfo <- function (file, RG, path = NULL, sep="\t", quote="\"") 
{
    if (!is.null(path)) 
        file <- file.path(path, file)
    Table <- read.table(file, sep = sep, header = T, quote=quote, as.is=TRUE, fill=TRUE)
    Chr <- as.character(Table$Chr)
    indX <- which(Chr == "X" | Chr == "x")
    indY <- which(Chr == "Y" | Chr == "y")
    Chr[indX] <- 23
    Chr[indY] <- 24
    Position <- Table$Position
    RG$genes <- data.frame(RG$genes, Position, Chr = as.numeric(Chr))
    RG
}


log2ratios <- function(x) {
  if(!is.null(x$M.observed)){
    matrix(as.matrix(x$M.observed), nrow = nrow(x$M.observed), ncol = ncol(x$M.observed), byrow = FALSE, dimnames = dimnames(x$M.observed))
  }
  else{
    matrix(as.matrix(x$M), nrow = nrow(x$M), ncol = ncol(x$M), byrow = FALSE, dimnames = dimnames(x$M))
  }
}


