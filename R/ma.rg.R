MA.RG <- function (object, bc.method = "subtract", offset = 0, design = c(1, ncol(object))) 
{
    if (is.null(object$R) || is.null(object$G)) 
        stop("Object doesn't contain R and G components")
    object <- backgroundCorrect(object, method = bc.method, offset = offset)
    R <- object$R
    G <- object$G
    R[R <= 0] <- NA
    G[G <= 0] <- NA
    R <- log(R, 2)
    G <- log(G, 2)
    object$R <- object$G <- object$Rb <- object$Gb <- object$other <- NULL
    object$M <- as.matrix((R - G))
##    if(is.null(design)){
##      design = c(rep(1, ncol(object$M)))
##    }
    for(i in 1:length(design)){
      temp <- design[i]* object$M[,i]
      object$M[,i] <- temp
    }
    object$A <- as.matrix((R + G)/2)
    new("MAList", unclass(object))
}

normalizeWithinArrays <- 
function (object, layout = object$printer, method = "printtiploess", 
    weights = object$weights, span = 0.3, iterations = 4, controlspots = NULL, 
    df = 5, robust = "M", bc.method = "subtract", offset = 0, design = c(1, ncol(object))) 
{
    if (!is(object, "MAList")) 
        object <- MA.RG(object, bc.method = bc.method, offset = offset, design = design)
    choices <- c("none", "median", "loess", "printtiploess", 
        "composite", "robustspline")
    method <- match.arg(method, choices)
    if (method == "none") 
        return(object)
    narrays <- ncol(object$M)
    if (method == "median") {
        for (j in 1:narrays) object$M[, j] <- object$M[, j] - 
            median(object$M[, j], na.rm = TRUE)
        return(object)
    }
    switch(method, loess = {
        for (j in 1:narrays) {
            y <- object$M[, j]
            x <- object$A[, j]
            w <- weights[, j]
            object$M[, j] <- loessFit(y, x, w, span = span, iterations = iterations)$residuals
        }
    }, printtiploess = {
        if (is.null(layout)) 
            stop("Layout argument not specified")
        ngr <- layout$ngrid.r
        ngc <- layout$ngrid.c
        nspots <- layout$nspot.r * layout$nspot.c
        for (j in 1:narrays) {
            spots <- 1:nspots
            for (gridr in 1:ngr) for (gridc in 1:ngc) {
                y <- object$M[spots, j]
                x <- object$A[spots, j]
                w <- weights[spots, j]
                object$M[spots, j] <- loessFit(y, x, w, span = span, 
                  iterations = iterations)$residuals
                spots <- spots + nspots
            }
        }
    }, composite = {
        if (is.null(layout)) 
            stop("Layout argument not specified")
        if (is.null(controlspots)) 
            stop("controlspots argument not specified")
        ntips <- layout$ngrid.r * layout$ngrid.c
        nspots <- layout$nspot.r * layout$nspot.c
        for (j in 1:narrays) {
            y <- object$M[, j]
            x <- object$A[, j]
            w <- weights[, j]
            f <- is.finite(y) & is.finite(x) & is.finite(w)
            y[!f] <- NA
            fit <- loess(y ~ x, weights = w, span = span, subset = controlspots, 
                na.action = na.exclude, degree = 0, surface = "direct", 
                family = "symmetric", trace.hat = "approximate", 
                iterations = iterations)
            alpha <- global <- y
            global[f] <- predict(fit, newdata = x[f])
            alpha[f] <- (rank(x[f]) - 1)/sum(f)
            spots <- 1:nspots
            for (tip in 1:ntips) {
                y <- object$M[spots, j]
                x <- object$A[spots, j]
                w <- weights[spots, j]
                local <- loessFit(y, x, w, span = span, iterations = iterations)$fitted
                object$M[spots, j] <- object$M[spots, j] - alpha[spots] * 
                  global[spots] - (1 - alpha[spots]) * local
                spots <- spots + nspots
            }
        }
    }, robustspline = {
        if (is.null(layout)) 
            stop("Layout argument not specified")
        for (j in 1:narrays) object$M[, j] <- normalizeRobustSpline(object$M[, 
            j], object$A[, j], layout, df = df, method = robust)
    })
    object
}
