"two.states" <-
function(data,covariates,iterlim,var.fixed=FALSE){

data <- as.matrix(data)

state <- 2

mu <- vector(length=state)
Sigma <- vector(length=state)
S <- vector(length=state)
emis.prob <- matrix(nrow=state,ncol=nrow(data),b=TRUE)
alpha <- matrix(nrow=state,ncol=nrow(data),b=TRUE)
alphahat <- matrix(nrow=state,ncol=nrow(data),b=TRUE)
denom <- vector(length=nrow(data))
prior <- vector(length=state)
dist <- vector()


fr.two.het <- function(x) {
  mu[1] <- x[1]
  mu[2] <- x[2]
  Sigma[1] <- x[3]
  Sigma[2] <- x[4]
  prior <- x[5]
  eta <- x[6]
  zeta <- x[7]
  omega <- x[8]
  S[1] <- exp(Sigma[1])
  if (var.fixed == FALSE) S[2] <- exp(Sigma[2]) else {S[2] <- exp(Sigma[1])}
  if (prior < 150) pr1 <- exp(prior)/(1+exp(prior)) else pr1 <- 1
  if (eta < 150) p1 <- exp(eta)/(1 + exp(eta)) else p1 <- 1
  if (zeta < 150) p2 <- exp(zeta)/(1 + exp(zeta)) else p2 <- 1
  if (omega < 50) rate1 <- exp(omega) else rate1 <- exp(50)
  gammaA <- matrix(ncol=2,nrow=2,b=T)
  gammaB <- matrix(ncol=2,nrow=2,b=T)
  gammaC <- matrix(ncol=2,nrow=2,b=T)
  gammaA <- matrix(c(1 - p1, p1, p2, 1 - p2 ), ncol = 2, b = TRUE)
  gammaB <- matrix(c(p1, -p1, -p2, p2), ncol = 2, b = TRUE)
  
  for (j in 1:nrow(data))
  {for (k in 1:2)
     { if (S[k] > 0.001) emis.prob[k,j] <-  dnorm(data[j,1], mean = mu[k], sd = S[k], log = FALSE) else {if (data[j,1] >= mu[k]*(1-0.001) & data[j,1] <= mu[k]*(1+0.001)) emis.prob[k,j] <- 1 else {emis.prob[k,j] <- 0}}}}
  alpha[1,1] <- pr1*emis.prob[1,1]
  alpha[2,1] <- (1 - pr1)*emis.prob[2,1]

  alphahat[1,1] <- alpha[1,1]/(sum(alpha[1,1],alpha[2,1]))
  alphahat[2,1] <- alpha[2,1]/(sum(alpha[1,1],alpha[2,1]))

  for (t in 2:nrow(data)){
    gammaC <- gammaA + exp(-((covariates[t-1,1]^rate1)*(prod(covariates[t-1,-1]))))*gammaB 
    for (i in 1:2)
           {alpha[i,t] <- (alphahat[1,t-1]*gammaC[1,i] +  alphahat[2,t-1]*gammaC[2,i]) * emis.prob[i,t]}
    for (i in 1:2){
            alphahat[i,t] <- alpha[i,t]/(sum(alpha[1,t],alpha[2,t]))
          }
  }

  for (t in 1:nrow(data)){
    denom[t] <- sum(alpha[1,t],alpha[2,t])
  }
  
 -(sum(log(denom)))
}




init.mean.two <- pam(data[,1],2)$medoids

init.var.two <- vector()

if (var.fixed == FALSE){
if (length(pam(data[,1],2)$data[pam(data[,1],2)$clustering==1]) > 1)  init.var.two[1] <- log(sqrt(var(pam(data[,1],2)$data[pam(data[,1],2)$clustering==1]))) else init.var.two[1] <- log(0.5)
if (length(pam(data[,1],2)$data[pam(data[,1],2)$clustering==2]) > 1)  init.var.two[2] <- log(sqrt(var(pam(data[,1],2)$data[pam(data[,1],2)$clustering==2]))) else init.var.two[2] <- log(0.5)} else {
init.var.two[1] <- log(sqrt(var(data[,1])))
init.var.two[2] <- log(sqrt(var(data[,1])))}

het.two <- nlm(fr.two.het ,c(init.mean.two[,1],init.var.two,-1,-3.6,-3.6,0),  iterlim=iterlim)

het.two
}

