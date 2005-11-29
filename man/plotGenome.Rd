\name{plotGenome}
\alias{plotGenome}
\alias{maxna}
\alias{minna}
\alias{sample.names}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plots the genome}
\description{
  Basic plot of the log2 ratios for each array ordered along the genome.
}
\usage{
plotGenome(input, array = 1, naut = 22, Y = TRUE, X
                 = TRUE, status, values, pch, col, cex, chrominfo =
                 chrominfo.basepair, yScale = c(-2, 2), samplenames =
                 (colnames(input))[array], ylb = "Log2Ratio",
                 title.font = 2, chrom.to.plot = 1, Z = FALSE, xlower =
                 0, xupper = clone.genomepos[sum(clone.genomepos >
                 0)]/1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{an object of class \code{\link[limma]{MAList}} or \code{\link[snapCGH:SegList]{SegList}}}
  \item{array}{integer of the array (sample) to be plotted.}
  \item{naut}{number of autosomes in the organism}
  \item{Y}{TRUE if chromosome Y is to be plotted, FALSE otherwise}
  \item{X}{TRUE if chromosome X is to be plotted, FALSE otherwise}
    \item{Z}{TRUE if you want to plot and individual chromosome (in
      which case \code{chrom.to.plot} specifies the chromosome to plot), FALSE otherwise}
%  \item{data}{a matrix containing values to use for plotting. defaults
%    to the \code{log2.ratios(input)}.}
  \item{status}{character vector giving the control status of each spot on the array, of same length as the number of rows of \code{log2.ratios(input)}. If omitted, all points are plotted in the default color, symbol and size.}
  \item{values}{character vector giving values of \code{status} to be highlighted on the plot. Defaults to unique values of \code{status}. Ignored if there is no \code{status} vector.}
  \item{pch}{vector or list of plotting characters. Default to integer code 16. Ignored is there is no \code{status} vector.}
  \item{col}{numeric or character vector of colors, of the same length as \code{values}. Defaults to \code{1:length(values)}. Ignored if there is no \code{status} vector.}
  \item{cex}{numeric vector of plot symbol expansions, of the the same
    length as \code{values}. Defaults to 0.2 for the most common
    \code{status} value and 1 for the others. Ignored if there is no
    \code{status} vector.}
  \item{chrominfo}{a chromosomal information associated with the mapping of the data.}
  \item{yScale}{Minimum y-scale to use for plotting. Scale is expanded
    if any of the values exceed the positive or negative limit.}
  \item{chrom.to.plot}{Vector of integers specifying chromosomes to plot
    if \code{Z=TRUE}}
  \item{samplenames}{sample names.}
  \item{ylb}{label for the Y-axis.}
  \item{title.font}{specify the font of the title; see
    \code{\link[base:cex]{cex}}}
  \item{xlower}{lower limit for the x-axis if \code{Z=TRUE}}
  \item{xupper}{upper limit for the x-axis if \code{Z=TRUE}}
}
\Details{
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
%plotGenome(MA)
}
\keyword{hplot}% at least one, from doc/KEYWORDS
