"find.param.one" <-
function(output.optim)
  {
    object <- list()
    object$mu <- output.optim$estimate[1]
    object$sigma <- output.optim$esitmate[2]
    object$minus.logLikelihood <- output.optim$minimum
    object$convergence <- output.optim$code
    object
  }

