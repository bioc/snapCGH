"find.param.five" <-
function(output.optim,var.fixed)
  {
    object <- list()
    object$mu <- c(output.optim$estimate[1],output.optim$estimate[2],output.optim$estimate[3],output.optim$estimate[4],output.optim$estimate[5])
    if (var.fixed == FALSE) {object$sigma <- c(exp(output.optim$estimate[6]),exp(output.optim$estimate[7]),exp(output.optim$estimate[8]),exp(output.optim$estimate[9]),exp(output.optim$estimate[10]))} else {
    object$sigma <- c(exp(output.optim$estimate[6]),exp(output.optim$estimate[6]),exp(output.optim$estimate[6]),exp(output.optim$estimate[6]),exp(output.optim$estimate[6]))}

    if (output.optim$estimate[11] < 150) pr1 <- exp(output.optim$estimate[11])/(1 + exp(output.optim$estimate[11])) else pr1 <- 1
    if (output.optim$estimate[12] < 150) pr2 <- (1 - pr1)*exp(output.optim$estimate[12])/(1 + exp(output.optim$estimate[12])) else pr2 <- (1 - pr1)
    if (output.optim$estimate[13] < 150) pr3 <- (1 - pr1 - pr2)*exp(output.optim$estimate[13])/(1 + exp(output.optim$estimate[13])) else pr3 <- (1 - pr1 - pr2)
    if (output.optim$estimate[14] < 150) pr4 <- (1 - pr1 - pr2 - pr3)*exp(output.optim$estimate[14])/(1 + exp(output.optim$estimate[14])) else pr4 <- (1 - pr1 - pr2 - pr3)
    
    if (output.optim$estimate[15] < 150) p1 <- exp(output.optim$estimate[15])/(1 + exp(output.optim$estimate[15])) else p1 <- 1
    if (output.optim$estimate[16] < 150) p2 <- (1 - p1)*exp(output.optim$estimate[16])/(1 + exp(output.optim$estimate[16])) else p2 <- (1 - p1)
    if (output.optim$estimate[17] < 150) p3 <- (1 - p1 - p2)*exp(output.optim$estimate[17])/(1 + exp(output.optim$estimate[17])) else p3 <- (1 - p1 - p2)
    if (output.optim$estimate[18] < 150) p4 <- (1-p1-p2-p3)*exp(output.optim$estimate[18])/(1 + exp(output.optim$estimate[18])) else p4 <- (1-p1-p2-p3)
    
    if (output.optim$estimate[19] < 150) p5 <- exp(output.optim$estimate[19])/(1 + exp(output.optim$estimate[19])) else p5 <- 1
    if (output.optim$estimate[20] < 150) p6 <- (1 - p5)*exp(output.optim$estimate[20])/(1 + exp(output.optim$estimate[20])) else p6 <- (1 - p5)
    if (output.optim$estimate[21] < 150) p7 <- (1 - p5 - p6)*exp(output.optim$estimate[21])/(1 + exp(output.optim$estimate[21])) else p7 <- (1 - p5 - p6)
    if (output.optim$estimate[22] < 150) p8 <- (1-p5-p6-p7)*exp(output.optim$estimate[22])/(1 + exp(output.optim$estimate[22])) else p8 <- (1-p5-p6-p7)
    
    if (output.optim$estimate[23] < 150) p9 <- exp(output.optim$estimate[23])/(1 + exp(output.optim$estimate[23])) else p9 <- 1
    if (output.optim$estimate[24] < 150) p10 <- (1 - p9)*exp(output.optim$estimate[24])/(1 + exp(output.optim$estimate[24])) else p10 <- (1 - p9)
    if (output.optim$estimate[25] < 150) p11 <- (1-p9- p10)*exp(output.optim$estimate[25])/(1 + exp(output.optim$estimate[25])) else p11 <- (1 - p9 - p10)
    if (output.optim$estimate[26] < 150) p12 <- (1-p9-p10-p11)*exp(output.optim$estimate[26])/(1 + exp(output.optim$estimate[26])) else p12 <- (1-p9-p10-p11)
    
    if (output.optim$estimate[27] < 150) p13 <- exp(output.optim$estimate[27])/(1 + exp(output.optim$estimate[27])) else p13 <- 1
    if (output.optim$estimate[28] < 150) p14 <- (1 - p13)*exp(output.optim$estimate[28])/(1 + exp(output.optim$estimate[28])) else p14 <- (1 - p13)
    if (output.optim$estimate[29] < 150) p15 <- (1 -p13-p14)*exp(output.optim$estimate[29])/(1 + exp(output.optim$estimate[29])) else p15 <- (1 -p13-p14)
    if (output.optim$estimate[30] < 150) p16 <- (1-p13-p14-p15)*exp(output.optim$estimate[30])/(1 + exp(output.optim$estimate[30])) else p16 <- (1-p13-p14-p15)
    
    if (output.optim$estimate[31] < 150) p17 <- exp(output.optim$estimate[31])/(1 + exp(output.optim$estimate[31])) else p17 <- 1
    if (output.optim$estimate[32] < 150) p18 <- (1 - p17)*exp(output.optim$estimate[32])/(1 + exp(output.optim$estimate[32])) else p18 <- (1 - p17)
    if (output.optim$estimate[33] < 150) p19 <- (1-p17-p18)*exp(output.optim$estimate[33])/(1 + exp(output.optim$estimate[33])) else p19 <- (1-p17-p18)
    if (output.optim$estimate[34] < 150) p20 <- (1-p17-p18-p19)*exp(output.optim$estimate[34])/(1 + exp(output.optim$estimate[34])) else p20 <- (1-p17-p18-p19)
    
    object$prior <- c(pr1,pr2,pr3,pr4,1-pr1-pr2-pr3-pr4)
    if (output.optim$estimate[35] < 50) object$rate1 <- exp(output.optim$estimate[35]) else object$rate1 <- exp(50)
    object$LH.trans <- matrix(c(1-p1-p2-p3-p4, p1, p2, p3, p4, p5, 1-p5-p6-p7-p8, p6, p7, p8, p9, p10, 1-p9-p10-p11-p12, p11, p12, p13, p14, p15, 1-p13-p14-p15-p16, p16, p17, p18, p19, p20, 1 - p17-p18-p19-p20),ncol=5,b=T)
    object$RH.trans <- matrix(c(p1+p2+p3+p4, -p1, -p2, -p3, -p4, -p5, p5+p6+p7+p8, -p6, -p7, -p8,-p9, -p10, p9+p10+p11+p12,-p11,-p12,-p13,-p14, -p15,p13+p14+p15+p16,-p16,-p17,-p18,-p19,-p20,p17+p18+p19+p20),ncol=5,b=T)
    object$minus.logLikelihood <- output.optim$minimum
    object$convergence <- output.optim$code
    object
  }

