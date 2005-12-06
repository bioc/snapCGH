"one.state" <-
function(data,iterlim){

data <- as.matrix(data)

state <- 1

mu <- vector(length=state)
Sigma <- vector(length=state)
S <- vector(length=state)
emis.prob <- matrix(nrow=state,ncol=nrow(data),b=TRUE)
alpha <- matrix(nrow=state,ncol=nrow(data),b=TRUE)
alphahat <- matrix(nrow=state,ncol=nrow(data),b=TRUE)
denom <- vector(length=nrow(data))
prior <- vector(length=state)
dist <- vector()

# The function for calculating the likelihood

fr.one <- function(x) {
  mu[1] <- x[1]
  Sigma[1] <- x[2]

  S[1] <- exp(Sigma[1])
  
  for (j in 1:nrow(data))
  {for (k in 1:1)
##K is number of states, j is the observation
     {emis.prob[k,j] <-  dnorm(data[j,1], mean = mu[k], sd = S[k], log = FALSE) }
 }
  alpha[1,1] <- emis.prob[1,1]

  alphahat[1,1] <- alpha[1,1]/(sum(alpha[1,1]))

  for (t in 2:nrow(data)){
    for (i in 1:1)
           {alpha[i,t] <- alphahat[1,t-1]* emis.prob[i,t]}
    for (i in 1:1)
      {alphahat[i,t] <- alpha[i,t]/(sum(alpha[1,t]))}
  }

  for (t in 1:nrow(data)){
    denom[t] <- sum(alpha[1,t])
  }
  
 -(sum(log(denom)))
 }


# We now optimise the function for each chromosome

init.mean <- mean(data)

homogeneous.one <- nlm(fr.one , c(init.mean,0.5), iterlim=iterlim)

homogeneous.one
}

