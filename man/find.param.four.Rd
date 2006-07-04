\name{find.param.four}
\alias{find.param.four}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Yields output when there are 4 underlying states}
\description{
  This function provides the output (means/variances, transition matrix, likelihood) when the heterogeneous HMM is fitted with four underlying states. 
  It is a workhorse of the process.data function. 
}
\usage{
find.param.four(output.optim,var.fixed)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{output.optim}{ The output from fitting a heterogeneous HMM when there are four underlying states }
  \item{var.fixed}{ Logical variable - TRUE if you want the variance to be tied across states. Defaults to FALSE }
}
\value{
  This function outputs the state means/variances, transition matrices, rate parameters, maximised likelihood and convergence information provided
  by fitting a heterogeneous HMM with four underlying states.}
\author{ John Marioni }


\keyword{ misc}
