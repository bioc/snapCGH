"three.states" <-
function(data,covariates,iterlim,var.fixed=FALSE){

data <- as.matrix(data)

state <- 3

mu <- vector(length=state)
Sigma <- vector(length=state)
S <- vector(length=state)
emis.prob <- matrix(nrow=state,ncol=nrow(data),b = T)
alpha <- matrix(nrow=state,ncol=nrow(data),b=TRUE)
alphahat <- matrix(nrow=state,ncol=nrow(data),b=TRUE)
denom <- vector(length=nrow(data))
prior <- vector(length=state)
dist <- vector()

fr.three.het <- function(x) {
  mu[1] <- x[1]
  mu[2] <- x[2]
  mu[3] <- x[3]
  Sigma[1] <- x[4]
  Sigma[2] <- x[5]
  Sigma[3] <- x[6]
  prior[1] <- x[7]
  prior[2] <- x[8]
  eta <- x[9]
  zeta <- x[10]
  theta <- x[11]
  beta <- x[12]
  gamma <- x[13]
  xi <- x[14]
  omega <- x[15]
  S[1] <- exp(Sigma[1])
  if (var.fixed == FALSE) {S[2] <- exp(Sigma[2])
                           S[3] <- exp(Sigma[3])} else {S[2] <- exp(Sigma[1])
                                                        S[3] <- exp(Sigma[1])}

  if (eta < 150) p1 <- exp(eta)/(1 + exp(eta)) else p1 <- 1
  if (zeta < 150) p2 <- (1 - p1)*(exp(zeta)/(1+exp(zeta))) else p2 <- (1-p1)
  if (theta < 150) p3 <- exp(theta)/(1 + exp(theta)) else p3 <- 1
  if (beta < 150) p4 <- (1 - p3)*(exp(beta)/(1+exp(beta))) else p4 <- (1 - p3)
  if (gamma < 150) p5 <- exp(gamma)/(1 + exp(gamma)) else p5 <- 1
  if (xi < 150) p6 <- (1 - p5)*(exp(xi)/(1+exp(xi))) else p6 <- (1-p5)
  if (prior[1] < 150) pr1 <- exp(prior[1])/(1 + exp(prior[1])) else pr1 <- 1
  if (prior[2] < 150) pr2 <- (1 - pr1)*(exp(prior[2]))/(1 + exp(prior[2])) else pr2 <- (1-pr1)
  if (omega) rate1 <- exp(omega) else rate1 <- exp(50)

  gammaA <- matrix(ncol=3,nrow=3,b=T)
  gammaB <- matrix(ncol=3,nrow=3,b=T)
  gammaC <- matrix(ncol=3,nrow=3,b=T)
  gammaA <- matrix(c(1 - p1 - p2, p1, p2, p3, 1 - p3 - p4,  p4, p5, p6, 1 - p5 - p6), ncol = 3, b = TRUE)
  gammaB <- matrix(c(p1 + p2, -p1, -p2, -p3, p3 + p4, -p4, -p5, -p6, p5 + p6), ncol = 3, b = TRUE)
  
  for (j in 1:nrow(data))
  {for (k in 1:3)
     { if (S[k] > 0.001) emis.prob[k,j] <-  dnorm(data[j,1], mean = mu[k], sd = S[k], log = FALSE) else {if (data[j,1] >= mu[k]*(1-0.01) & data[j,1] <= mu[k]*(1+0.01)) emis.prob[k,j] <- 1 else {emis.prob[k,j] <- 0}}}}
  
  alpha[1,1] <- pr1*emis.prob[1,1]
  alpha[2,1] <- pr2*emis.prob[2,1]
  alpha[3,1] <- (1 - pr1 - pr2)*emis.prob[3,1]

  alphahat[1,1] <- alpha[1,1]/(sum(alpha[1,1],alpha[2,1],alpha[3,1]))
  alphahat[2,1] <- alpha[2,1]/(sum(alpha[1,1],alpha[2,1],alpha[3,1]))
  alphahat[3,1] <- alpha[3,1]/(sum(alpha[1,1],alpha[2,1],alpha[3,1]))

  for (t in 2:nrow(data)){
    gammaC <- gammaA + exp(-((covariates[t-1,1]^rate1)*(prod(covariates[t-1,-1]))))*gammaB 
    for (i in 1:3)
           {alpha[i,t] <- (alphahat[1,t-1]*gammaC[1,i] +  alphahat[2,t-1]*gammaC[2,i] + alphahat[3,t-1]*gammaC[3,i]) * emis.prob[i,t]}
    for (i in 1:3)
      {alphahat[i,t] <- alpha[i,t]/(sum(alpha[1,t],alpha[2,t],alpha[3,t]))
     }
  }

  for (t in 1:nrow(data))
    {denom[t] <- sum(alpha[1,t],alpha[2,t],alpha[3,t])
   }
 -(sum(log(denom)))
 }

# Initialisation

init.mean.three <- vector()
init.mean.three <- pam(data[,1],3)$medoids

init.var.three <- vector()


if (var.fixed == FALSE){
if (length(pam(data[,1],3)$data[pam(data[,1],3)$clustering==1]) > 1)  init.var.three[1] <- log(sqrt(var(pam(data[,1],3)$data[pam(data[,1],3)$clustering==1]))) else init.var.three[1] <- log(0.5)
if (length(pam(data[,1],3)$data[pam(data[,1],3)$clustering==2]) > 1)  init.var.three[2] <- log(sqrt(var(pam(data[,1],3)$data[pam(data[,1],3)$clustering==2]))) else init.var.three[2] <- log(0.5)
if (length(pam(data[,1],3)$data[pam(data[,1],3)$clustering==3]) > 1)  init.var.three[3] <- log(sqrt(var(pam(data[,1],3)$data[pam(data[,1],3)$clustering==3]))) else init.var.three[3] <- log(0.5)} else {
init.var.three[1] <- log(sqrt(var(data[,1])))
init.var.three[2] <- log(sqrt(var(data[,1])))
init.var.three[3] <- log(sqrt(var(data[,1])))}


# Optimisation

het.three <- nlm(fr.three.het ,c(init.mean.three[,1],init.var.three,-0.7,-0.7,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,0),iterlim=iterlim)
  
het.three
}

