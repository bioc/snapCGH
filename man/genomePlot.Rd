\name{genomePlot}
\alias{genomePlot}
\alias{sample.names}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plots the genome}
\description{
  Basic plot of the log2 ratios for each array ordered along the genome.
}
\usage{
genomePlot(input, array = 1, naut = 22, Y = FALSE, X
                 = FALSE, main = NA, status, values, pch, cex, col, chrominfo =
                 chrominfo.Mb, ylim = c(-2, 2), ylb = "Log2Ratio",
                 chrom.to.plot = NA, xlim = c(0,NA), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{an object of class \code{\link[limma]{MAList}} or \code{\link[snapCGH:SegList]{SegList}}}
  \item{array}{integer of the array (sample) to be plotted.}
  \item{naut}{number of autosomes in the organism}
  \item{Y}{TRUE if chromosome Y is to be plotted, FALSE otherwise}
  \item{X}{TRUE if chromosome X is to be plotted, FALSE otherwise}
  \item{main}{Provides the title of the plot}	
  \item{status}{character vector giving the control status of each spot on the array, of same length as the number of rows of \code{log2ratios(input)}. If omitted, all points are plotted in the default color, symbol and size.}
  \item{values}{character vector giving values of \code{status} to be highlighted on the plot. Defaults to unique values of \code{status}. Ignored if there is no \code{status} vector.}
  \item{pch}{vector or list of plotting characters. Default to integer code 16. Ignored is there is no \code{status} vector.}
  \item{col}{numeric or character vector of colors, of the same length as \code{values}. Defaults to \code{1:length(values)}. Ignored if there is no \code{status} vector.}
  \item{cex}{numeric vector of plot symbol expansions, of the the same
    length as \code{values}. Defaults to 0.2 for the most common
    \code{status} value and 1 for the others. Ignored if there is no
    \code{status} vector.}
  \item{chrominfo}{a chromosomal information associated with the mapping of the data.}
  \item{ylim}{Minimum y-scale to use for plotting.}
  \item{chrom.to.plot}{Specify which chromosome to plot}
  \item{ylb}{label for the Y-axis.}
  \item{xlim}{limits for the x-axis}
  \item{...}{Any other parameters}
}
\details{
  The status vector is intended to specify the control status of each
  spot, for example "gene", "ratio control", "house keeping gene",
  "buffer" and so on. The vector is usually computed using the function
  controlStatus and a spot-types file. However the function may be used
  to highlight any subset of spots.
}
\seealso{
  \code{\link[limma]{MAList}}
  \code{\link[snapCGH:SegList]{SegList}}
}
\author{John Marioni}
\examples{

%#reading the spot types file
%spotTypes <- readSpotTypes("SpotTypes.txt")
%MA$genes$Status <- controlStatus(spotTypes, MA)
%genomePlot(MA)
}
\keyword{hplot}% at least one, from doc/KEYWORDS
