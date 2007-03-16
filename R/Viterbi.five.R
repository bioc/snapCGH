"Viterbi.five" <- function(y,BFGS.output,BFGS.trans.mat)  {
  
  K <- 5
  n <- length(y)

  psi <- matrix(nrow=K,ncol=n,b=TRUE)
  phi <- matrix(nrow=K,ncol=n,b=TRUE)
  q.star <- vector()

  for (i in 1:K){
    phi[i,1] <- max(-100000,log(BFGS.output$prior[i]),na.rm=TRUE) + log(dnorm(y[1],BFGS.output$mu[i], BFGS.output$sigma[i], log = FALSE))
  }

  for (i in 1:K){
    psi[i,1] <- 0
  }

  for (t in 2:n) {
    for (k in 1:K){
      phi[k,t] <- max(phi[1,t-1] + max(-100000,log(BFGS.trans.mat[[t-1]][1,k]),na.rm=TRUE), phi[2,t-1] + max(-100000,log(BFGS.trans.mat[[t-1]][2,k]),na.rm=TRUE),phi[3,t-1] + max(-100000,log(BFGS.trans.mat[[t-1]][3,k]),na.rm=TRUE),phi[4,t-1] + max(-100000,log(BFGS.trans.mat[[t-1]][4,k]),na.rm=TRUE),phi[5,t-1] + max(-100000,log(BFGS.trans.mat[[t-1]][5,k]),na.rm=TRUE)) + log(dnorm(y[t],BFGS.output$mu[k], BFGS.output$sigma[k], log = FALSE))
    }
  }
  


for (t in 2:n) {
  for (k in 1:K){
    if ((phi[1,t-1] + max(-100000,log(BFGS.trans.mat[[t-1]][1,k]),na.rm=TRUE)) > max(phi[2,t-1] + max(-100000,log(BFGS.trans.mat[[t-1]][2,k]),na.rm=TRUE),phi[3,t-1] + max(-100000,log(BFGS.trans.mat[[t-1]][3,k]),na.rm=TRUE),phi[4,t-1] + max(-100000,log(BFGS.trans.mat[[t-1]][4,k]),na.rm=TRUE),phi[5,t-1] + max(-100000,log(BFGS.trans.mat[[t-1]][5,k]),na.rm=TRUE))) psi[k,t] <- 1 else {
      if ((phi[2,t-1] + max(-100000,log(BFGS.trans.mat[[t-1]][2,k]),na.rm=TRUE)) > max(phi[3,t-1] + max(-100000,log(BFGS.trans.mat[[t-1]][3,k]),na.rm=TRUE),phi[4,t-1] + max(-100000,log(BFGS.trans.mat[[t-1]][4,k]),na.rm=TRUE),phi[5,t-1] + max(-100000,log(BFGS.trans.mat[[t-1]][5,k]),na.rm=TRUE))) psi[k,t] <- 2 else {
        if ((phi[3,t-1] + max(-100000,log(BFGS.trans.mat[[t-1]][3,k]),na.rm=TRUE)) > (max(phi[4,t-1] + max(-100000,log(BFGS.trans.mat[[t-1]][4,k]),na.rm=TRUE),phi[5,t-1] + max(-100000,log(BFGS.trans.mat[[t-1]][5,k]),na.rm=TRUE)))) psi[k,t] <- 3 else {
          if ((phi[4,t-1] + max(-100000,log(BFGS.trans.mat[[t-1]][4,k]),na.rm=TRUE)) > (phi[5,t-1] + max(-100000,log(BFGS.trans.mat[[t-1]][5,k]),na.rm=TRUE)))   psi[k,t] <- 4 else {psi[k,t] <- 5}}}}
  }
}



  P.star <- max(phi[1,n],phi[2,n],phi[3,n],phi[4,n],phi[5,n])
  if (phi[1,n] > max(phi[2,n],phi[3,n],phi[4,n],phi[5,n])) q.star[n] <- 1 else {if (phi[2,n] > max(phi[3,n],phi[4,n],phi[5,n])) q.star[n] <- 2 else {if (phi[3,n] > max(phi[4,n],phi[5,n])) q.star[n] <- 3 else {if (phi[4,n] > phi[5,n]) q.star[n] <- 4 else {q.star[n] <- 5}}}}
  
  for (t in ((n-1):1)){
    q.star[t] <- psi[q.star[t+1],t+1]
  }
  
  q.star
}

