"find.param.five" <-
function(output.optim,var.fixed)
  {
    object <- list()
    object$mu <- c(output.optim[[1]][1],output.optim[[1]][2],output.optim[[1]][3],output.optim[[1]][4],output.optim[[1]][5])
    if (var.fixed == FALSE) {object$sigma <- c(exp(output.optim[[1]][6]),exp(output.optim[[1]][7]),exp(output.optim[[1]][8]),exp(output.optim[[1]][9]),exp(output.optim[[1]][10]))} else {
    object$sigma <- c(exp(output.optim[[1]][6]),exp(output.optim[[1]][6]),exp(output.optim[[1]][6]),exp(output.optim[[1]][6]),exp(output.optim[[1]][6]))}

    if (output.optim[[1]][11] < 150) pr1 <- exp(output.optim[[1]][11])/(1 + exp(output.optim[[1]][11])) else pr1 <- 1
    if (output.optim[[1]][12] < 150) pr2 <- (1 - pr1)*exp(output.optim[[1]][12])/(1 + exp(output.optim[[1]][12])) else pr2 <- (1 - pr1)
    if (output.optim[[1]][13] < 150) pr3 <- (1 - pr1 - pr2)*exp(output.optim[[1]][13])/(1 + exp(output.optim[[1]][13])) else pr3 <- (1 - pr1 - pr2)
    if (output.optim[[1]][14] < 150) pr4 <- (1 - pr1 - pr2 - pr3)*exp(output.optim[[1]][14])/(1 + exp(output.optim[[1]][14])) else pr4 <- (1 - pr1 - pr2 - pr3)
    
    if (output.optim[[1]][15] < 150) p1 <- exp(output.optim[[1]][15])/(1 + exp(output.optim[[1]][15])) else p1 <- 1
    if (output.optim[[1]][16] < 150) p2 <- (1 - p1)*exp(output.optim[[1]][16])/(1 + exp(output.optim[[1]][16])) else p2 <- (1 - p1)
    if (output.optim[[1]][17] < 150) p3 <- (1 - p1 - p2)*exp(output.optim[[1]][17])/(1 + exp(output.optim[[1]][17])) else p3 <- (1 - p1 - p2)
    if (output.optim[[1]][18] < 150) p4 <- (1-p1-p2-p3)*exp(output.optim[[1]][18])/(1 + exp(output.optim[[1]][18])) else p4 <- (1-p1-p2-p3)
    
    if (output.optim[[1]][19] < 150) p5 <- exp(output.optim[[1]][19])/(1 + exp(output.optim[[1]][19])) else p5 <- 1
    if (output.optim[[1]][20] < 150) p6 <- (1 - p5)*exp(output.optim[[1]][20])/(1 + exp(output.optim[[1]][20])) else p6 <- (1 - p5)
    if (output.optim[[1]][21] < 150) p7 <- (1 - p5 - p6)*exp(output.optim[[1]][21])/(1 + exp(output.optim[[1]][21])) else p7 <- (1 - p5 - p6)
    if (output.optim[[1]][22] < 150) p8 <- (1-p5-p6-p7)*exp(output.optim[[1]][22])/(1 + exp(output.optim[[1]][22])) else p8 <- (1-p5-p6-p7)
    
    if (output.optim[[1]][23] < 150) p9 <- exp(output.optim[[1]][23])/(1 + exp(output.optim[[1]][23])) else p9 <- 1
    if (output.optim[[1]][24] < 150) p10 <- (1 - p9)*exp(output.optim[[1]][24])/(1 + exp(output.optim[[1]][24])) else p10 <- (1 - p9)
    if (output.optim[[1]][25] < 150) p11 <- (1-p9- p10)*exp(output.optim[[1]][25])/(1 + exp(output.optim[[1]][25])) else p11 <- (1 - p9 - p10)
    if (output.optim[[1]][26] < 150) p12 <- (1-p9-p10-p11)*exp(output.optim[[1]][26])/(1 + exp(output.optim[[1]][26])) else p12 <- (1-p9-p10-p11)
    
    if (output.optim[[1]][27] < 150) p13 <- exp(output.optim[[1]][27])/(1 + exp(output.optim[[1]][27])) else p13 <- 1
    if (output.optim[[1]][28] < 150) p14 <- (1 - p13)*exp(output.optim[[1]][28])/(1 + exp(output.optim[[1]][28])) else p14 <- (1 - p13)
    if (output.optim[[1]][29] < 150) p15 <- (1 -p13-p14)*exp(output.optim[[1]][29])/(1 + exp(output.optim[[1]][29])) else p15 <- (1 -p13-p14)
    if (output.optim[[1]][30] < 150) p16 <- (1-p13-p14-p15)*exp(output.optim[[1]][30])/(1 + exp(output.optim[[1]][30])) else p16 <- (1-p13-p14-p15)
    
    if (output.optim[[1]][31] < 150) p17 <- exp(output.optim[[1]][31])/(1 + exp(output.optim[[1]][31])) else p17 <- 1
    if (output.optim[[1]][32] < 150) p18 <- (1 - p17)*exp(output.optim[[1]][32])/(1 + exp(output.optim[[1]][32])) else p18 <- (1 - p17)
    if (output.optim[[1]][33] < 150) p19 <- (1-p17-p18)*exp(output.optim[[1]][33])/(1 + exp(output.optim[[1]][33])) else p19 <- (1-p17-p18)
    if (output.optim[[1]][34] < 150) p20 <- (1-p17-p18-p19)*exp(output.optim[[1]][34])/(1 + exp(output.optim[[1]][34])) else p20 <- (1-p17-p18-p19)
    
    object$prior <- c(pr1,pr2,pr3,pr4,1-pr1-pr2-pr3-pr4)
    if (output.optim[[1]][35] < 50) object$rate1 <- exp(output.optim[[1]][35]) else object$rate1 <- exp(50)
    object$LH.trans <- matrix(c(1-p1-p2-p3-p4, p1, p2, p3, p4, p5, 1-p5-p6-p7-p8, p6, p7, p8, p9, p10, 1-p9-p10-p11-p12, p11, p12, p13, p14, p15, 1-p13-p14-p15-p16, p16, p17, p18, p19, p20, 1 - p17-p18-p19-p20),ncol=5,b=TRUE)
    object$RH.trans <- matrix(c(p1+p2+p3+p4, -p1, -p2, -p3, -p4, -p5, p5+p6+p7+p8, -p6, -p7, -p8,-p9, -p10, p9+p10+p11+p12,-p11,-p12,-p13,-p14, -p15,p13+p14+p15+p16,-p16,-p17,-p18,-p19,-p20,p17+p18+p19+p20),ncol=5,b=TRUE)
    object$minus.logLikelihood <- output.optim$val
    # object$convergence <- output.optim$code
    object
  }

