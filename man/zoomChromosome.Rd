\name{zoomChromosome}
\alias{zoomChromosome}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Interactive plot of a single chromsome}
\description{
  Plot splitting the screen into two.  The top windows displays the
  entire chromosome, whilst the bottom plots a selected region.  The plot
  is interactive allowing the user to click twice on a chromosome in the upper
  plot and have it the region between the two clicks displayed below.
}
\usage{
zoomChromosome(..., array = 1, chrom.to.plot, colors = NULL, chrominfo =
chrominfo.Mb, ylim = c(-2, 2))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{...}{Objects of type \code{\link[limma:malist]{MAList}} or
	  \code{\link[snapCGH:SegList]{SegList}}}
	\item{array}{Specify which array should be plotted.}
	\item{chrom.to.plot}{Which chromosome should be plotted}
	\item{colors}{Vector specify the colors for each of the SegLists}
	\item{chrominfo}{chromosomal information associated with the
	  mapping of the data.}
	\item{ylim}{Specify the limits of the y-axis}
}
\details{
  If colors is unspecified then all SegLists passed to this function will
  be plotted in blue.  Since this makes it quite hard to tell which is
  which it is highly recommended to specify the colors vector if more
  than one object is being passe to this function.
}

\author{Mike Smith}

\keyword{hplot}
