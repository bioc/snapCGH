\name{find.param.two}
\alias{find.param.two}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Yields output when there are 2 underlying states}
\description{
  This function provides the output (means/variances, transition matrix, likelihood) when the heterogeneous HMM is fitted with two underlying states. 
  It is a workhorse of the process.data function. 
}
\usage{
find.param.two(output.optim,var.fixed)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{output.optim}{ The output from fitting a heterogeneous HMM when there are two underlying states }
  \item{var.fixed}{ Logical variable - TRUE if you want to tie the variance across states. Defaults to FALSE }
}
\value{
 This function outputs the state means/variances, transition matrices, rate parameters, maximised likelihood and convergence information provided
  by fitting a heterogeneous HMM with two underlying states.}
\author{ John Marioni }

\keyword{ misc}
