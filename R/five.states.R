"five.states" <-
function(data,covariates,iterlim,var.fixed=FALSE){

data <- as.matrix(data)

state <- 5

mu <- vector(length=state)
Sigma <- vector(length=state)
S <- vector(length=state)
emis.prob <- matrix(nrow=state,ncol=nrow(data),b=T)
alpha <- matrix(nrow=state,ncol=nrow(data),b=TRUE)
alphahat <- matrix(nrow=state,ncol=nrow(data),b=TRUE)
denom <- vector(length=nrow(data))
prior <- vector(length=state)
dist <- vector()

fr.five.het <- function(x) {
  mu[1] <- x[1]
  mu[2] <- x[2]
  mu[3] <- x[3]
  mu[4] <- x[4]
  mu[5] <- x[5]
  Sigma[1] <- x[6]
  Sigma[2] <- x[7]
  Sigma[3] <- x[8]
  Sigma[4] <- x[9]
  Sigma[5] <- x[10]
  prior[1] <- x[11]
  prior[2] <- x[12]
  prior[3] <- x[13]
  prior[4] <- x[14]
  eta <- x[15]
  zeta <- x[16]
  nu <- x[17]
  omikron <- x[18] 
  theta <- x[19]
  beta <- x[20]
  phi <- x[21]
  kappa <- x[22]
  gamma <- x[23]
  delt <- x[24]
  epsilon <- x[25]
  tau <- x[26]
  lambda <- x[27]
  rho <- x[28]
  xi <- x[29]
  iota <- x[30]
  chi <- x[31]
  upsilon <- x[32]
  psi <- x[33]
  aux <- x[34]
  omega <- x[35]
  S[1] <- exp(Sigma[1])
  if (var.fixed == FALSE){
  S[2] <- exp(Sigma[2])
  S[3] <- exp(Sigma[3])
  S[4] <- exp(Sigma[4])
  S[5] <- exp(Sigma[5])} else {
  S[2] <- exp(Sigma[1])
  S[3] <- exp(Sigma[1])
  S[4] <- exp(Sigma[1])
  S[5] <- exp(Sigma[1])}
  if (eta < 150) p1 <- exp(eta)/(1 + exp(eta)) else p1 <- 1
  if (zeta < 150) p2 <- (1 - p1)*(exp(zeta)/(1+exp(zeta))) else p2 <- (1-p1)
  if (nu < 150) p3 <- (1 - p1 - p2)*(exp(nu)/(1+exp(nu))) else p3 <- (1 - p1 - p2)
  if (omikron < 150) p4 <- (1 - p1 - p2 - p3)*(exp(omikron)/(1+exp(omikron))) else p4 <- (1 - p1 - p2 - p3)
  if (theta < 150) p5 <- exp(theta)/(1 + exp(theta)) else p5 <- 1
  if (beta < 150) p6 <- (1 - p5)*(exp(beta)/(1+exp(beta))) else p6 <- (1-p5)
  if (phi < 150) p7 <- (1 - p5 - p6)*(exp(phi)/(1+exp(phi))) else p7 <- (1 - p5 - p6)
  if (kappa < 150) p8 <- (1 - p5 - p6 - p7)*(exp(kappa)/(1+exp(kappa))) else p8 <- (1 - p5 - p6 - p7)
  if (gamma < 150) p9 <- exp(gamma)/(1 + exp(gamma)) else p9 <- 1
  if (delt < 150) p10 <- (1 - p9)*(exp(delt)/(1+exp(delt))) else p10 <- (1-p9)
  if (epsilon < 150) p11 <- (1 - p9 - p10)*(exp(epsilon)/(1+exp(epsilon))) else p11 <- (1 - p9 - p10)
  if (tau < 150) p12 <- (1 - p9 - p10 - p11)*(exp(tau)/(1+exp(tau))) else p12 <- (1 - p9 - p10 - p11)
  if (lambda < 150) p13 <- exp(lambda)/(1 + exp(lambda)) else p13 <- 1
  if (rho < 150) p14 <- (1 - p13)*(exp(rho)/(1+exp(rho))) else p14 <- (1 - p13)
  if (xi < 150) p15 <- (1 - p13 - p14)*(exp(xi)/(1+exp(xi))) else p15 <- (1 - p13 - p14)
  if (iota < 150) p16 <- (1 - p13 - p14 - p15)*(exp(iota)/(1+exp(iota))) else p16 <- (1 - p13 - p14 - p15)
  if (chi < 150) p17 <- exp(chi)/(1 + exp(chi)) else p17 <- 1
  if (upsilon < 150) p18 <- (1 - p17)*(exp(upsilon)/(1+exp(upsilon))) else p18 <- (1-p17)
  if (psi < 150) p19 <- (1 - p17 - p18)*(exp(psi)/(1+exp(psi))) else p19 <- (1 - p17 - p18)
  if (aux < 150) p20 <- (1 - p17 - p18 - p19)*(exp(aux)/(1+exp(aux))) else p20 <- (1 - p17 - p18 - p19)
  if (prior[1] < 150) pr1 <- exp(prior[1])/(1 + exp(prior[1])) else pr1 <- 1
  if (prior[2] < 150) pr2 <- (1 - pr1)*(exp(prior[2]))/(1 + exp(prior[2])) else pr2 <- (1-pr1)
  if (prior[3] < 150) pr3 <- (1 - pr1 - pr2)*(exp(prior[3])/(1 + exp(prior[3]))) else pr3 <- (1 - pr1 - pr2)
  if (prior[4] < 150) pr4 <- (1 - pr1 - pr2 - pr3)*(exp(prior[4])/(1 + exp(prior[4]))) else pr4 <- (1 - pr1 - pr2 - pr3)
  pr5 <- 1 - pr1 - pr2 - pr3 - pr4
  if (omega < 50) rate1 <- exp(omega) else rate1 <- exp(50)
  
  gammaA <- matrix(ncol=5,nrow=5,b=T)
  gammaB <- matrix(ncol=5,nrow=5,b=T)
  gammaC <- matrix(ncol=5,nrow=5,b=T)
  gammaA <- matrix(c(abs(1 - p1 - p2 - p3 - p4), p1, p2, p3, p4, p5, abs(1 - p5 - p6 - p7 - p8),  p6, p7, p8, p9, p10, abs(1 - p9 -p10 - p11 -p12), p11, p12, p13, p14, p15, abs(1 - p13 - p14 - p15 - p16), p16, p17, p18, p19, p20, abs(1 - p17 - p18 - p19 - p20)), ncol = 5, b = TRUE)
  gammaB <- matrix(c(p1 + p2 + p3 + p4, -p1, -p2, -p3, -p4, -p5, p5 + p6 + p7 + p8, -p6, -p7, -p8, -p9, -p10, p9 + p10 + p11 + p12, -p11, -p12, -p13, -p14, -p15, p13 + p14 + p15 + p16, -p16, -p17, -p18, -p19, -p20, p17 + p18 + p19 + p20), ncol = 5, b = TRUE)
  
  for (j in 1:nrow(data))
  {for (k in 1:5)
     { if (S[k] > 0.001) emis.prob[k,j] <-  dnorm(data[j,1], mean = mu[k], sd = S[k], log = FALSE) else {if (data[j,1] >= mu[k]*(1-0.001) & data[j,1] <= mu[k]*(1+0.001)) emis.prob[k,j] <- 1 else {emis.prob[k,j] <- 0}}}}
  
  alpha[1,1] <- pr1*emis.prob[1,1]
  alpha[2,1] <- pr2*emis.prob[2,1]
  alpha[3,1] <- pr3*emis.prob[3,1]
  alpha[4,1] <- pr4*emis.prob[4,1]
  alpha[5,1] <- pr5*emis.prob[5,1]

  alphahat[1,1] <- alpha[1,1]/(sum(alpha[1,1],alpha[2,1],alpha[3,1],alpha[4,1],alpha[5,1]))
  alphahat[2,1] <- alpha[2,1]/(sum(alpha[1,1],alpha[2,1],alpha[3,1],alpha[4,1],alpha[5,1]))
  alphahat[3,1] <- alpha[3,1]/(sum(alpha[1,1],alpha[2,1],alpha[3,1],alpha[4,1],alpha[5,1]))
  alphahat[4,1] <- alpha[4,1]/(sum(alpha[1,1],alpha[2,1],alpha[3,1],alpha[4,1],alpha[5,1]))
  alphahat[5,1] <- alpha[5,1]/(sum(alpha[1,1],alpha[2,1],alpha[3,1],alpha[4,1],alpha[5,1]))
  
  for (t in 2:nrow(data)){
     gammaC <- gammaA + exp(-((covariates[t-1,1]^rate1)*(prod(covariates[t-1,-1]))))*gammaB 
    for (j in 1:5)
      {alpha[j,t] <- (alphahat[1,t-1]*gammaC[1,j] +  alphahat[2,t-1]*gammaC[2,j] + alphahat[3,t-1]*gammaC[3,j] + alphahat[4,t-1]*gammaC[4,j] + alphahat[5,t-1]*gammaC[5,j]) * emis.prob[j,t]}
    for (j in 1:5){
      alphahat[j,t] <- alpha[j,t]/(sum(alpha[1,t],alpha[2,t],alpha[3,t],alpha[4,t],alpha[5,t]))
    }
  }

  for (t in 1:nrow(data)){
    denom[t] <- sum(alpha[1,t],alpha[2,t],alpha[3,t],alpha[4,t],alpha[5,t])
  }
  -(sum(log(denom)))
 }


# Initialisation

init.mean.five <- vector()
init.mean.five <- pam(data[,1],5)$medoids

init.var.five <- list()

init.var.five <- vector()

if (var.fixed == FALSE){
if (length(pam(data[,1],5)$data[pam(data[,1],5)$clustering==1]) > 1)  init.var.five[1] <- log(sqrt(var(pam(data[,1],5)$data[pam(data[,1],5)$clustering==1]))) else init.var.five[1] <- log(0.5)
if (length(pam(data[,1],5)$data[pam(data[,1],5)$clustering==2]) > 1)  init.var.five[2] <- log(sqrt(var(pam(data[,1],5)$data[pam(data[,1],5)$clustering==2]))) else init.var.five[2] <- log(0.5)
if (length(pam(data[,1],5)$data[pam(data[,1],5)$clustering==3]) > 1)  init.var.five[3] <- log(sqrt(var(pam(data[,1],5)$data[pam(data[,1],5)$clustering==3]))) else init.var.five[3] <- log(0.5)
if (length(pam(data[,1],5)$data[pam(data[,1],5)$clustering==4]) > 1)  init.var.five[4] <- log(sqrt(var(pam(data[,1],5)$data[pam(data[,1],5)$clustering==4]))) else init.var.five[4] <- log(0.5)
if (length(pam(data[,1],5)$data[pam(data[,1],5)$clustering==5]) > 1)  init.var.five[5] <- log(sqrt(var(pam(data[,1],5)$data[pam(data[,1],5)$clustering==5]))) else init.var.five[5] <- log(0.5)} else {
init.var.five[1] <- log(sqrt(var(data[,1])))
init.var.five[2] <- log(sqrt(var(data[,1])))
init.var.five[3] <- log(sqrt(var(data[,1])))
init.var.five[4] <- log(sqrt(var(data[,1])))
init.var.five[5] <- log(sqrt(var(data[,1])))}


# Optimisation

het.five <- nlm(fr.five.het, c(init.mean.five[,1],init.var.five,-0.7,-0.7,-0.7,-0.7,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,0), iterlim=iterlim)

het.five
}

