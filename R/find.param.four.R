"find.param.four" <-
function(output.optim,var.fixed)
  {
    object <- list()
    object$mu <- c(output.optim$x[1],output.optim$x[2],output.optim$x[3],output.optim$x[4])
    if (var.fixed == FALSE) {object$sigma <- c(exp(output.optim$x[5]),exp(output.optim$x[6]),exp(output.optim$x[7]),exp(output.optim$x[8]))} else {
    object$sigma <- c(exp(output.optim$x[5]),exp(output.optim$x[5]),exp(output.optim$x[5]),exp(output.optim$x[5]))}

    if (output.optim$x[9] < 150) pr1 <- exp(output.optim$x[9])/(1 + exp(output.optim$x[9])) else pr1 <- 1
    if (output.optim$x[10] < 150) pr2 <- (1 - pr1)*exp(output.optim$x[10])/(1 + exp(output.optim$x[10])) else pr2 <- (1-pr1)
    if (output.optim$x[11] < 150) pr3 <- (1-pr1-pr2)*exp(output.optim$x[11])/(1 + exp(output.optim$x[11])) else pr3 <- (1-pr1-pr2)
    
    if (output.optim$x[12] < 150) p1 <- exp(output.optim$x[12])/(1 + exp(output.optim$x[12])) else p1 <- 1
    if (output.optim$x[13] < 150) p2 <- (1 - p1)*exp(output.optim$x[13])/(1 + exp(output.optim$x[13])) else p2 <- (1 - p1)
    if (output.optim$x[14] < 150) p3 <- (1 - p1 - p2)*exp(output.optim$x[14])/(1 + exp(output.optim$x[14])) else p3 <- (1 - p1 - p2)
    
    if (output.optim$x[15] < 150) p4 <- exp(output.optim$x[15])/(1 + exp(output.optim$x[15])) else p4 <- 1
    if (output.optim$x[16] < 150) p5 <- (1 - p4)*exp(output.optim$x[16])/(1 + exp(output.optim$x[16])) else p5 <- (1 - p4)
    if (output.optim$x[17] < 150) p6 <- (1 - p4 - p5)*exp(output.optim$x[17])/(1 + exp(output.optim$x[17])) else p6 <- (1 - p4 - p5)
    
    if (output.optim$x[18] < 150) p7 <- exp(output.optim$x[18])/(1 + exp(output.optim$x[18])) else p7 <- 1
    if (output.optim$x[19] < 150) p8 <- (1 - p7)*exp(output.optim$x[19])/(1 + exp(output.optim$x[19])) else p8 <- (1 - p7)
    if (output.optim$x[20] < 150) p9 <- (1 - p7 - p8)*exp(output.optim$x[20])/(1 + exp(output.optim$x[20])) else p9 <- (1 - p7 - p8)
    
    if (output.optim$x[21] < 150) p10 <- exp(output.optim$x[21])/(1 + exp(output.optim$x[21])) else p10 <- 1
    if (output.optim$x[22] < 150) p11 <- (1 - p10)*exp(output.optim$x[22])/(1 + exp(output.optim$x[22])) else p11 <- (1 - p10)
    if (output.optim$x[23] < 150) p12 <- (1 - p10-p11)*exp(output.optim$x[23])/(1 + exp(output.optim$x[23])) else p12 <- (1-p10-p11)
    
    object$prior <- c(pr1,pr2,pr3,1-pr1-pr2-pr3)
    if (output.optim$x[24] < 700) object$rate1 <- exp(output.optim$x[24]) else object$rate1 <- exp(50)
    object$LH.trans <- matrix(c(abs(1-p1-p2-p3), p1, p2, p3, p4, abs(1-p4-p5-p6), p5, p6, p7, p8, abs(1-p7-p8-p9), p9, p10, p11, p12, abs(1-p10-p11-p12)),ncol=4,b=TRUE)
    object$RH.trans <- matrix(c(p1+p2+p3, -p1, -p2, -p3, -p4, p4 + p5 + p6, -p5, -p6, -p7,-p8,p7+p8+p9,-p9,-p10,-p11,-p12,p10+p11+p12),ncol=4,b=TRUE)
    object$minus.logLikelihood <- output.optim$val
    # object$convergence <- output.optim$code
    object
  }

