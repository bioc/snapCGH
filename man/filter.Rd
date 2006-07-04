\name{filterClones}
\alias{filterClones}
\title{Filter clones from sample}
\description{Function for filtering clones via a user defined function.}
\usage{
	filterClones(MA, filterFunc, ...)
}
\arguments{
	\item{MA}{An object of class \code{\link[limma:malist]{MAList}}}
	\item{filterFunc}{A user specified function that accepts an MAList 
	and returns the indices of the clones to be removed.}
	\item{...}{Additional arguments to be passed to the filter function.}
}
\details{
	Any clones identified by the filter function are turned into NA's.
	These are then removed or imputed within the \code{\link{processCGH}} function.
}

\author{Mike Smith}

\keyword{methods}
