\name{find.param.five}
\alias{find.param.five}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Yields the output in a model with five underlying states}
\description{
  This function is a workhorse of the process.data function. It outputs state means/variances and transitions matrices in the model with five states.
}
\usage{
find.param.five(output.optim,var.fixed)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{output.optim}{ output of the optimisation with 5 underlying states}
  \item{var.fixed}{Logical variable - TRUE if you want to fix the variance to be the same across states}
}
\value{
  Outputs the mean/variance, transition matrix, maximised likelihood and convergence information}
\author{ John Marioni }

\keyword{ misc}
