\name{Viterbi.three}
\alias{Viterbi.three}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ A scaled Viterbi algorithm for allocating clones to one of two underlying states.}
\description{
  A work horse of the fit.model function. It uses a scaled version of the Viterbi algorithm to allocate clones to one of three underlying states as fitted 
  using a heterogeneous HMM.
}
\usage{
Viterbi.three(y, BFGS.output, BFGS.trans)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ the data to be allocated to states}
  \item{BFGS.output}{ The output obtained from the find.param.three function }
  \item{BFGS.trans}{  A list of the heterogeneous transition matrices }
}
\value{
  A vector of numbers indicating to which state clones are allocated to.
}
\author{ John Marioni}

\keyword{ misc}
