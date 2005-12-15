\name{removeByWeights}
\alias{removeByWeights}
\title{Remove clones based on a weights matrix}
\description{An example function to be used by the filterClones method.
	This function takes an MA list, a weights matrix and a threshold 
	and returns the indices of any clones with weight below the 
	threshold.
}
\usage{
	removeByWeights(MA, weights=MA$weights, threshold = 0.2)
}
\arguments{
	\item{MA}{An object of class \code{\link[limma:malist]{MAList}}}
	\item{weights}{A matrix with the same dimensions as MA containing
		weight information.}
	\item{threshold}{Threshold value.  Any clones with weight below this 
		are removed.}
}

\author{Mike Smith}

\seealso{
	\code{\link{filterClones}}
}

\keyword{methods}