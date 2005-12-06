\name{one.state}
\alias{one.state}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function for finding the log-likelihood and maximising it when there is 1 underlying state}
\description{
  This function calculates the log-likelihood using a scaled procdeure and the forward-backward equations when there is 1 underlying state. It then 
  determines initial parameters before using the nlm function to optimise the likelihood. It is a work horse of the fit.model funtion.
}
\usage{
one.state(data, iterlim)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ The data to be segmented}
  \item{iterlim}{ The maximum number of iterations to be performed}
}
\value{
  The same format as returned by the nlm function
}
\author{ John Marioni}


\keyword{models}

