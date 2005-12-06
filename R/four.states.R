"four.states" <-
function(data,covariates,iterlim,var.fixed=FALSE){

data <- as.matrix(data)

state <- 4

mu <- vector(length=state)
Sigma <- vector(length=state)
S <- vector(length=state)
emis.prob <- matrix(nrow=state,ncol=nrow(data),b = T)
alpha <- matrix(nrow=state,ncol=nrow(data),b=TRUE)
alphahat <- matrix(nrow=state,ncol=nrow(data),b=TRUE)
denom <- vector(length=nrow(data))
prior <- vector(length=state)
dist <- vector()

fr.four.het <- function(x) {
  mu[1] <- x[1]
  mu[2] <- x[2]
  mu[3] <- x[3]
  mu[4] <- x[4]
  Sigma[1] <- x[5]
  Sigma[2] <- x[6]
  Sigma[3] <- x[7]
  Sigma[4] <- x[8]
  prior[1] <- x[9]
  prior[2] <- x[10]
  prior[3] <- x[11]
  eta <- x[12]
  zeta <- x[13]
  nu <- x[14]
  theta <- x[15]
  beta <- x[16]
  phi <- x[17]
  gamma <- x[18]
  delt <- x[19]
  epsilon <- x[20]
  lambda <- x[21]
  rho <- x[22]
  xi <- x[23]
  omega <- x[24]
  S[1] <- exp(Sigma[1])
  if (var.fixed == FALSE) {S[2] <- exp(Sigma[2])
                           S[3] <- exp(Sigma[3])
                           S[4] <- exp(Sigma[4])} else {
                           S[2] <- exp(Sigma[1])
                           S[3] <- exp(Sigma[1])
                           S[4] <- exp(Sigma[1])}
  if (eta < 150) p1 <- exp(eta)/(1 + exp(eta)) else p1 <- 1
  if (zeta < 150) p2 <- (1 - p1)*(exp(zeta)/(1+exp(zeta))) else p2 <- 1 - p1
  if (nu < 150) p3 <- (1 - p1 - p2)*(exp(nu)/(1+exp(nu))) else p3 <-  (1-p1-p2)
  if (theta < 150) p4 <- exp(theta)/(1 + exp(theta)) else p4 <- 1
  if (beta < 150) p5 <- (1 - p4)*(exp(beta)/(1+exp(beta))) else p5 <- 1 - p4
  if (phi < 150) p6 <- (1 - p4 - p5)*(exp(phi)/(1+exp(phi))) else p6 <- (1-p4-p5)
  if (gamma < 150) p7 <- exp(gamma)/(1 + exp(gamma)) else p7 <- 1
  if (delt < 150) p8 <- (1 - p7)*(exp(delt)/(1+exp(delt))) else p8 <- 1 - p7
  if (epsilon < 150) p9 <- (1 - p7 - p8)*(exp(epsilon)/(1+exp(epsilon))) else p9 <- (1-p7-p8)
  if (lambda < 150) p10 <- exp(lambda)/(1+exp(lambda)) else p10 <- 1
  if (rho < 150) p11 <- (1 - p10)*(exp(rho)/(1+exp(rho))) else p11 <- 1 - p10
  if (xi < 150) p12 <- (1 - p10 - p11)*(exp(xi)/(1+exp(xi))) else p12 <- (1-p10-p11)
  if (prior[1] < 150) pr1 <- exp(prior[1])/(1 + exp(prior[1])) else pr1 <- 1
  if (prior[2] < 150) pr2 <- (1 - pr1)*(exp(prior[2]))/(1 + exp(prior[2])) else pr2 <- 1 - pr1
  if (prior[3] < 150) pr3 <- (1 - pr1 - pr2)*(exp(prior[3])/(1 + exp(prior[3]))) else pr3 <- (1-pr1-pr2)
  rate1 <- exp(omega)
  
  gammaA <- matrix(ncol=4,nrow=4,b=T)
  gammaB <- matrix(ncol=4,nrow=4,b=T)
  gammaC <- matrix(ncol=4,nrow=4,b=T)
  gammaA <- matrix(c(1 - p1 - p2 - p3, p1, p2, p3, p4, 1 - p4 - p5 - p6,  p5, p6, p7, p8, 1 - p7 - p8 - p9, p9, p10, p11, p12, 1 - p10 - p11 - p12), ncol = 4, b = TRUE)
  gammaB <- matrix(c(p1 + p2 + p3, -p1, -p2, -p3, -p4, p4 + p5 + p6, -p5, -p6, -p7, -p8, p7 + p8 + p9, -p9, -p10, -p11, -p12, p10 + p11 + p12), ncol = 4, b = TRUE)
  
  for (j in 1:nrow(data))
  {for (k in 1:4)
     { if (S[k] > 0.001) emis.prob[k,j] <-  dnorm(data[j,1], mean = mu[k], sd = S[k], log = FALSE) else {if (data[j,1] >= (mu[k]*(1-0.001)) & data[j,1] <= mu[k]*(1+0.001)) emis.prob[k,j] <- 1 else {emis.prob[k,j] <- 0}}}}
     
  alpha[1,1] <- pr1*emis.prob[1,1]
  alpha[2,1] <- pr2*emis.prob[2,1]
  alpha[3,1] <- pr3*emis.prob[3,1]
  alpha[4,1] <- (1-pr1-pr2-pr3)*emis.prob[4,1]

  alphahat[1,1] <- alpha[1,1]/(sum(alpha[1,1],alpha[2,1],alpha[3,1],alpha[4,1]))
  alphahat[2,1] <- alpha[2,1]/(sum(alpha[1,1],alpha[2,1],alpha[3,1],alpha[4,1]))
  alphahat[3,1] <- alpha[3,1]/(sum(alpha[1,1],alpha[2,1],alpha[3,1],alpha[4,1]))
  alphahat[4,1] <- alpha[4,1]/(sum(alpha[1,1],alpha[2,1],alpha[3,1],alpha[4,1]))
  

  for (t in 2:nrow(data)){
    gammaC <- gammaA + exp(-((covariates[t-1,1]^rate1)*(prod(covariates[t-1,-1]))))*gammaB 
         for (i in 1:4)
           {alpha[i,t] <- (alphahat[1,t-1]*gammaC[1,i] +  alphahat[2,t-1]*gammaC[2,i] + alphahat[3,t-1]*gammaC[3,i] + alphahat[4,t-1]*gammaC[4,i]) * emis.prob[i,t]}
    for (i in 1:4)
      {alphahat[i,t] <- alpha[i,t]/(sum(alpha[1,t],alpha[2,t],alpha[3,t],alpha[4,t]))}
  }

  for (t in 1:nrow(data)){
    denom[t] <- sum(alpha[1,t],alpha[2,t],alpha[3,t],alpha[4,t])
  }

 -(sum(log(denom)))
 }

# Initialisation

init.mean.four <- vector()
init.mean.four <- pam(data[,1],4)$medoids

init.var.four <- list()

init.var.four <- vector()

if (var.fixed == FALSE){
if (length(pam(data[,1],4)$data[pam(data[,1],4)$clustering==1]) > 1)  init.var.four[1] <- log(sqrt(var(pam(data[,1],4)$data[pam(data[,1],4)$clustering==1]))) else init.var.four[1] <- log(0.5)
if (length(pam(data[,1],4)$data[pam(data[,1],4)$clustering==2]) > 1)  init.var.four[2] <- log(sqrt(var(pam(data[,1],4)$data[pam(data[,1],4)$clustering==2]))) else init.var.four[2] <- log(0.5)
if (length(pam(data[,1],4)$data[pam(data[,1],4)$clustering==3]) > 1)  init.var.four[3] <- log(sqrt(var(pam(data[,1],4)$data[pam(data[,1],4)$clustering==3]))) else init.var.four[3] <- log(0.5)
if (length(pam(data[,1],4)$data[pam(data[,1],4)$clustering==4]) > 1)  init.var.four[4] <- log(sqrt(var(pam(data[,1],4)$data[pam(data[,1],4)$clustering==4]))) else init.var.four[4] <- log(0.5)} else {
init.var.four[1] <- log(sqrt(var(data[,1])))
init.var.four[2] <- log(sqrt(var(data[,1])))
init.var.four[3] <- log(sqrt(var(data[,1])))
init.var.four[4] <- log(sqrt(var(data[,1])))}


# Optimisation

het.four <- nlm(fr.four.het, c(init.mean.four[,1],init.var.four,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,0),iterlim=iterlim)

het.four
}

