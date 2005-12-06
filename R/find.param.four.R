"find.param.four" <-
function(output.optim,var.fixed)
  {
    object <- list()
    object$mu <- c(output.optim$estimate[1],output.optim$estimate[2],output.optim$estimate[3],output.optim$estimate[4])
    if (var.fixed == FALSE) {object$sigma <- c(exp(output.optim$estimate[5]),exp(output.optim$estimate[6]),exp(output.optim$estimate[7]),exp(output.optim$estimate[8]))} else {
    object$sigma <- c(exp(output.optim$estimate[5]),exp(output.optim$estimate[5]),exp(output.optim$estimate[5]),exp(output.optim$estimate[5]))}

    if (output.optim$estimate[9] < 150) pr1 <- exp(output.optim$estimate[9])/(1 + exp(output.optim$estimate[9])) else pr1 <- 1
    if (output.optim$estimate[10] < 150) pr2 <- (1 - pr1)*exp(output.optim$estimate[10])/(1 + exp(output.optim$estimate[10])) else pr2 <- (1-pr1)
    if (output.optim$estimate[11] < 150) pr3 <- (1-pr1-pr2)*exp(output.optim$estimate[11])/(1 + exp(output.optim$estimate[11])) else pr3 <- (1-pr1-pr2)
    
    if (output.optim$estimate[12] < 150) p1 <- exp(output.optim$estimate[12])/(1 + exp(output.optim$estimate[12])) else p1 <- 1
    if (output.optim$estimate[13] < 150) p2 <- (1 - p1)*exp(output.optim$estimate[13])/(1 + exp(output.optim$estimate[13])) else p2 <- (1 - p1)
    if (output.optim$estimate[14] < 150) p3 <- (1 - p1 - p2)*exp(output.optim$estimate[14])/(1 + exp(output.optim$estimate[14])) else p3 <- (1 - p1 - p2)
    
    if (output.optim$estimate[15] < 150) p4 <- exp(output.optim$estimate[15])/(1 + exp(output.optim$estimate[15])) else p4 <- 1
    if (output.optim$estimate[16] < 150) p5 <- (1 - p4)*exp(output.optim$estimate[16])/(1 + exp(output.optim$estimate[16])) else p5 <- (1 - p4)
    if (output.optim$estimate[17] < 150) p6 <- (1 - p4 - p5)*exp(output.optim$estimate[17])/(1 + exp(output.optim$estimate[17])) else p6 <- (1 - p4 - p5)
    
    if (output.optim$estimate[18] < 150) p7 <- exp(output.optim$estimate[18])/(1 + exp(output.optim$estimate[18])) else p7 <- 1
    if (output.optim$estimate[19] < 150) p8 <- (1 - p7)*exp(output.optim$estimate[19])/(1 + exp(output.optim$estimate[19])) else p8 <- (1 - p7)
    if (output.optim$estimate[20] < 150) p9 <- (1 - p7 - p8)*exp(output.optim$estimate[20])/(1 + exp(output.optim$estimate[20])) else p9 <- (1 - p7 - p8)
    
    if (output.optim$estimate[21] < 150) p10 <- exp(output.optim$estimate[21])/(1 + exp(output.optim$estimate[21])) else p10 <- 1
    if (output.optim$estimate[22] < 150) p11 <- (1 - p10)*exp(output.optim$estimate[22])/(1 + exp(output.optim$estimate[22])) else p11 <- (1 - p10)
    if (output.optim$estimate[23] < 150) p12 <- (1 - p10-p11)*exp(output.optim$estimate[23])/(1 + exp(output.optim$estimate[23])) else p12 <- (1-p10-p11)
    
    object$prior <- c(pr1,pr2,pr3,1-pr1-pr2-pr3)
    if (output.optim$estimate[24] < 700) object$rate1 <- exp(output.optim$estimate[24]) else object$rate1 <- exp(50)
    object$LH.trans <- matrix(c(abs(1-p1-p2-p3), p1, p2, p3, p4, abs(1-p4-p5-p6), p5, p6, p7, p8, abs(1-p7-p8-p9), p9, p10, p11, p12, abs(1-p10-p11-p12)),ncol=4,b=T)
    object$RH.trans <- matrix(c(p1+p2+p3, -p1, -p2, -p3, -p4, p4 + p5 + p6, -p5, -p6, -p7,-p8,p7+p8+p9,-p9,-p10,-p11,-p12,p10+p11+p12),ncol=4,b=T)
    object$minus.logLikelihood <- output.optim$minimum
    object$convergence <- output.optim$code
    object
  }

