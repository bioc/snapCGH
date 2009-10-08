\name{zoomGenome}
\alias{zoomGenome}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Interactive plot of the whole genome}
\description{
  Plot splitting the screen into two.  The top windows displays the
  entire genome, whilst the bottom plots a single chromosome.  The plot
  is interactive allowing the user to click on a chromosome in the upper
  plot and have it displayed below.  Clicking to either side of the plot borders
  ends the interactivity.
}
\usage{
zoomGenome(..., array = 1, colors = NULL, chrominfo = chrominfo.Mb)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{...}{Objects of type \code{\link[snapCGH:SegList]{SegList}}}
	\item{array}{Specify which array should be plotted.}
	\item{colors}{Vector specify the colors for each of the SegLists}
	\item{chrominfo}{chromosomal information associated with the mapping of the data.}
}
\details{
  If colors is unspecified then all objects passed to this function will
  be plotted in blue.  Since this makes it quite hard to tell which is
  which it is highly recommended to specify the colors vector if more
  than one object is being passe to this function.
}

\author{Mike Smith}

\keyword{hplot}
