\name{plotSegmentationStates}
\alias{plotSegmentationStates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plotting the estimated hmm states and log2 ratios for each sample.
}
\description{
  This function displays the estimated hmm states and log2 ratios for
  each sample.
}
\usage{

plotSegmentationStates(segList, geList, array=1, chr = 1:length(unique(segList$genes$Chr)), maxChrom = 22,
                 chrominfo = chrominfo.basepair, yScale = c(-2, 2),
                 samplenames = sample.names(segList), xlower = 0,
                 xupper = chrominfo$length[chr[j]]/1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{segList}{object of class \code{\link[snapCGH:SegList]{SegList}}}
  \item{geList}{object of class \code{\link[snapCGH:GEList]{GEList}}}
  \item{array}{index of the array to be plotted}
  \item{chr}{vector of chromosomes to be plotted}
  \item{yScale}{specified scale for Y-exis}
  \item{maxChrom}{highest chromosome to show}
  \item{chrominfo}{a chromosomal information associated with the mapping of the data}
  \item{samplenames}{vector of sample names}
  \item{xlower}{integer value specifying the lower limit of the x-axis}
  \item{xupper}{integer value specifying the upper limit of the x-axis}
}
\details{
  Each chromosome is plotted on a separate  page and contains two
  figures. The top figure shows the observed log2ratios and the bottom
  figure shows predicted values for all clones but outliers which show
  observed values. The genomic events are indicated on both figures as
  following. The first clone after transition is indicated with solid
  blue line and the last clone after transitions is shown with dotted
  green line. Focal aberrations clones are colored orange,
  amplifications are colored red and outliers are yellow.
}
\references{Application of Hidden Markov Models to the analysis of the
  array CGH data, Fridlyand et.al., \emph{JMVA}, 2004
}
\author{Jane Fridlyand}
\seealso{
  \code{\link[snapCGH:SegList]{SegList}}
  \code{\link{fitHMM}}
  \code{\link{plotGenome}}
}
\examples{

%data(colorectal)
%plotSegmentationStates(colorectal, 1)

}
\keyword{hplot}% at least one, from doc/KEYWORDS
