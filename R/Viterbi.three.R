"Viterbi.three" <-
function(y,BFGS.output,BFGS.trans)
  {

K <- 3
n <- length(y)
  
psi <- matrix(nrow=K,ncol=n,byrow=TRUE)
phi <- matrix(nrow=K,ncol=n,byrow=TRUE)
q.star <- vector()

for (i in 1:K)
  {  phi[i,1] <- max(-100000,log(BFGS.output$prior[i]),na.rm=TRUE) + log(dnorm(y[1],BFGS.output$mu[i], BFGS.output$sigma[i], log = FALSE))
   }

for (i in 1:K)
  { psi[i,1] <- 0
  }

for (t in 2:n) {
  for (k in 1:K){
    phi[k,t] <- max(phi[1,t-1] + max(-100000,log(BFGS.trans[[t-1]][1,k]),na.rm=TRUE), phi[2,t-1] + max(-100000,log(BFGS.trans[[t-1]][2,k]),na.rm=TRUE),phi[3,t-1] + max(-100000,log(BFGS.trans[[t-1]][3,k]),na.rm=TRUE)) + log(dnorm(y[t],BFGS.output$mu[k], BFGS.output$sigma[k], log = FALSE))
  }
}



for (t in 2:n) {
  for (k in 1:K){
    if ((phi[1,t-1] + max(-100000,log(BFGS.trans[[t-1]][1,k]),na.rm=TRUE)) > (max(phi[2,t-1] + max(-100000,log(BFGS.trans[[t-1]][2,k]),na.rm=TRUE),phi[3,t-1] + max(-100000,log(BFGS.trans[[t-1]][3,k]),na.rm=TRUE)))) psi[k,t] <- 1 else {if ((phi[2,t-1] + max(-100000,log(BFGS.trans[[t-1]][2,k]),na.rm=TRUE)) > (phi[3,t-1] + max(-100000,log(BFGS.trans[[t-1]][3,k]),na.rm=TRUE))) psi[k,t] <- 2 else psi[k,t] <- 3}
  }
}


P.star <- max(phi[1,n],phi[2,n],phi[3,n])
  if (phi[1,n] > max(phi[2,n],phi[3,n])) q.star[n] <- 1 else {if (phi[2,n] > phi[3,n]) q.star[n] <- 2 else q.star[n] <- 3}

for (t in ((n-1):1)){
  q.star[t] <- psi[q.star[t+1],t+1]
}

q.star
}

