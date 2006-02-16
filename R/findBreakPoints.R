findBreakPoints <- function(states){
  bpoints <- c(1)
  for(i in 2:length(states)){
    if(states[i] != states[i-1])
      bpoints <- c(bpoints,i-1,i)
  }
  bpoints <- c(bpoints, length(states))
  bpoints
}


