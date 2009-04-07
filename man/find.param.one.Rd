\name{find.param.one}
\alias{find.param.one}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Yields output when there is 1 underlying states}
\description{
  This function provides the output (means/variances, transition matrix, likelihood) when the heterogeneous HMM is fitted with one underlying state. 
  It is a workhorse of the process.data function. 

}
\usage{
find.param.one(output.optim)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{output.optim}{ The output from fitting a heterogeneous HMM when there is one underlying states }
}
\value{
  This function outputs the state means/variances, transition matrices, rate parameters, maximised likelihood and convergence information provided
  by fitting a heterogeneous HMM with one underlying state.
}
\author{ John Marioni }

\keyword{ misc}

