\name{four.states}
\alias{four.states}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function for finding the log-likelihood and maximising it when there are 4 underlying states}
\description{
  This function calculates the log-likelihood using a scaled procdeure and the forward-backward equations when there are 4 underlying states. It then 
  determines initial parameters before using the nlm function to optimise the likelihood. It is a work horse of the fit.model funtion.
}
\usage{
four.states(data, covariates, iterlim,var.fixed=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ The data to be segmented}
  \item{covariates}{ The covariate information to be used in the segmentation }
  \item{iterlim}{  The maximum number of iterations to be performed }
  \item{var.fixed}{ Logical variable - TRUE if you want to tie the variance to be the same across states. Defaults to FALSE.}
}
\value{
  The same format as returned by the nlm function}
\author{ John Marioni }

\keyword{models}