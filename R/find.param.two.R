"find.param.two" <-
function(output.optim,var.fixed)
  {
    object <- list()
    object$mu <- c(output.optim$x[1],output.optim$x[2])
    if (var.fixed == FALSE) {object$sigma <- c(exp(output.optim$x[3]),exp(output.optim$x[4]))} else {
    object$sigma <- c(exp(output.optim$x[3]),exp(output.optim$x[3]))}
    if (output.optim$x[5] < 150) pr1 <- exp(output.optim$x[5])/(1 + exp(output.optim$x[5])) else pr1 <- 1
    if (output.optim$x[6] < 150) p1 <- exp(output.optim$x[6])/(1 + exp(output.optim$x[6])) else p1 <- 1
    if (output.optim$x[7] < 150) p2 <- exp(output.optim$x[7])/(1 + exp(output.optim$x[7])) else p2 <- 1
    object$prior <- c(pr1,1-pr1)
    if (output.optim$x[8] < 50) object$rate1 <- exp(output.optim$x[8]) else object$rate1 <- exp(50)
    object$LH.trans <- matrix(c(1 - p1, p1, p2, 1- p2),ncol=2,b=T)
    object$RH.trans <- matrix(c(p1, -p1, -p2, p2),ncol=2,b=T)
    object$minus.logLikelihood <- output.optim$val
    # object$convergence <- output.optim$code
    object
  }

