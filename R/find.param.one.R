"find.param.one" <-
function(output.optim)
  {
    object <- list()
    object$mu <- output.optim$x[1]
    object$sigma <- output.optim$x[2]
    object$minus.logLikelihood <- output.optim$val
    #object$convergence <- output.optim$code
    object
  }

