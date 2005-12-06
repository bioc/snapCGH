"find.param.two" <-
function(output.optim,var.fixed)
  {
    object <- list()
    object$mu <- c(output.optim$estimate[1],output.optim$estimate[2])
    if (var.fixed == FALSE) {object$sigma <- c(exp(output.optim$estimate[3]),exp(output.optim$estimate[4]))} else {
    object$sigma <- c(exp(output.optim$estimate[3]),exp(output.optim$estimate[3]))}
    if (output.optim$estimate[5] < 150) pr1 <- exp(output.optim$estimate[5])/(1 + exp(output.optim$estimate[5])) else pr1 <- 1
    if (output.optim$estimate[6] < 150) p1 <- exp(output.optim$estimate[6])/(1 + exp(output.optim$estimate[6])) else p1 <- 1
    if (output.optim$estimate[7] < 150) p2 <- exp(output.optim$estimate[7])/(1 + exp(output.optim$estimate[7])) else p2 <- 1
    object$prior <- c(pr1,1-pr1)
    if (output.optim$estimate[8] < 50) object$rate1 <- exp(output.optim$estimate[8]) else object$rate1 <- exp(50)
    object$LH.trans <- matrix(c(1 - p1, p1, p2, 1- p2),ncol=2,b=T)
    object$RH.trans <- matrix(c(p1, -p1, -p2, p2),ncol=2,b=T)
    object$minus.logLikelihood <- output.optim$minimum
    object$convergence <- output.optim$code
    object
  }

