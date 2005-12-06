"find.param.three" <-
function(output.optim,var.fixed)
  {
    object <- list()
    object$mu <- c(output.optim$estimate[1],output.optim$estimate[2],output.optim$estimate[3])
    if (var.fixed == FALSE) {object$sigma <- c(exp(output.optim$estimate[4]),exp(output.optim$estimate[5]),exp(output.optim$estimate[6]))} else {
    object$sigma <- c(exp(output.optim$estimate[4]),exp(output.optim$estimate[4]),exp(output.optim$estimate[4]))}
    if (output.optim$estimate[7] < 150) pr1 <- exp(output.optim$estimate[7])/(1 + exp(output.optim$estimate[7])) else pr1 <- 1
    if (output.optim$estimate[8] < 150) pr2 <- (1 - pr1)*exp(output.optim$estimate[8])/(1 + exp(output.optim$estimate[8])) else pr2 <- (1-pr1)
    if (output.optim$estimate[9] < 150) p1 <- exp(output.optim$estimate[9])/(1 + exp(output.optim$estimate[9])) else p1 <- 1
    if (output.optim$estimate[10] < 150) p2 <- (1 - p1)*exp(output.optim$estimate[10])/(1 + exp(output.optim$estimate[10])) else p2 <- (1-p1)
    if (output.optim$estimate[11] < 150) p3 <- exp(output.optim$estimate[11])/(1 + exp(output.optim$estimate[11])) else p3 <- 1
    if (output.optim$estimate[12] < 150) p4 <- (1 - p3)*exp(output.optim$estimate[12])/(1 + exp(output.optim$estimate[12])) else p4 <- (1-p3)
    if (output.optim$estimate[13] < 150) p5 <- exp(output.optim$estimate[13])/(1 + exp(output.optim$estimate[13])) else p5 <- 1
    if (output.optim$estimate[14] < 150) p6 <- (1 - p5)*exp(output.optim$estimate[14])/(1 + exp(output.optim$estimate[14])) else p6 <- (1-p5)
    object$prior <- c(pr1,pr2,1-pr1-pr2)
    if (output.optim$estimate[15] < 150) object$rate1 <- exp(output.optim$estimate[15]) else object$rate1 <- exp(50)
    object$LH.trans <- matrix(c(abs(1 - p1 - p2), p1, p2, p3, abs(1 - p3 - p4), p4, p5, p6, abs(1- p5 - p6)),ncol=3,b=T)
    object$RH.trans <- matrix(c(p1+p2, -p1, -p2, -p3, p3+p4, -p4, -p5, -p6, p5+p6),ncol=3,b=T)
    object$minus.logLikelihood <- output.optim$minimum
    object$convergence <- output.optim$code
    object
  }

