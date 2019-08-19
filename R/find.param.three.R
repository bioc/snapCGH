"find.param.three" <-
function(output.optim,var.fixed)
  {
    object <- list()
    object$mu <- c(output.optim$x[1],output.optim$x[2],output.optim$x[3])
    if (var.fixed == FALSE) {object$sigma <- c(exp(output.optim$x[4]),exp(output.optim$x[5]),exp(output.optim$x[6]))} else {
    object$sigma <- c(exp(output.optim$x[4]),exp(output.optim$x[4]),exp(output.optim$x[4]))}
    if (output.optim$x[7] < 150) pr1 <- exp(output.optim$x[7])/(1 + exp(output.optim$x[7])) else pr1 <- 1
    if (output.optim$x[8] < 150) pr2 <- (1 - pr1)*exp(output.optim$x[8])/(1 + exp(output.optim$x[8])) else pr2 <- (1-pr1)
    if (output.optim$x[9] < 150) p1 <- exp(output.optim$x[9])/(1 + exp(output.optim$x[9])) else p1 <- 1
    if (output.optim$x[10] < 150) p2 <- (1 - p1)*exp(output.optim$x[10])/(1 + exp(output.optim$x[10])) else p2 <- (1-p1)
    if (output.optim$x[11] < 150) p3 <- exp(output.optim$x[11])/(1 + exp(output.optim$x[11])) else p3 <- 1
    if (output.optim$x[12] < 150) p4 <- (1 - p3)*exp(output.optim$x[12])/(1 + exp(output.optim$x[12])) else p4 <- (1-p3)
    if (output.optim$x[13] < 150) p5 <- exp(output.optim$x[13])/(1 + exp(output.optim$x[13])) else p5 <- 1
    if (output.optim$x[14] < 150) p6 <- (1 - p5)*exp(output.optim$x[14])/(1 + exp(output.optim$x[14])) else p6 <- (1-p5)
    object$prior <- c(pr1,pr2,1-pr1-pr2)
    if (output.optim$x[15] < 150) object$rate1 <- exp(output.optim$x[15]) else object$rate1 <- exp(50)
    object$LH.trans <- matrix(c(abs(1 - p1 - p2), p1, p2, p3, abs(1 - p3 - p4), p4, p5, p6, abs(1- p5 - p6)),ncol=3,byrow=TRUE)
    object$RH.trans <- matrix(c(p1+p2, -p1, -p2, -p3, p3+p4, -p4, -p5, -p6, p5+p6),ncol=3,byrow=TRUE)
    object$minus.logLikelihood <- output.optim$val
    # object$convergence <- output.optim$code
    object
  }

