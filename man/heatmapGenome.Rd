\name{heatmapGenome}
\alias{heatmapGenome}
\alias{plotvalChrom.func}
\alias{plotValGenome}
\alias{plotValChrom}
\alias{plotChrom}
\alias{maPalette}
\alias{floor.func}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{clustering and heatmap}
\description{
  This function clusters samples and shows their heatmap
}
\usage{
heatmapGenome(input, response = as.factor(rep("All", ncol(input))),
                 chrominfo = chrominfo.Mb, cutoff = 1, lowCol =
                 "blue", highCol = "yellow", midCol = "white", ncolors =
                 50, byclass = FALSE, showaber = FALSE, amplif = 1,
                 homdel = -0.75, samplenames = colnames(input),
                 vecchrom = 1:22, titles = "Image Plot", methodS =
                 "ward", categoricalPheno = TRUE, CENTROMERE = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{object of class \code{\link[limma:malist]{MAList}} or
    \code{\link[snapCGH:SegList]{SegList}}}
  \item{response}{phenotype of interest. defaults to the same phenotype
    assigned to all samples}
  \item{chrominfo}{a chromosomal information associated with the mapping
    of the data}
  \item{cutoff}{maximum absolute value. all the values are floored to
    +/-cutoff depending on whether they are positive of
    negative. defaults to 1}
  \item{ncolors}{number of colors in the grid. input to
    \code{\link{maPalette}}. defaults to 50}
  \item{lowCol}{color for the low (negative) values. input to
    \code{\link{maPalette}}. defaults to "red"}
  \item{highCol}{color for the high (positive) values. input to
    \code{\link{maPalette}}. defaults to "green"}
  \item{midCol}{color for the values close to 0. input to
    \code{\link{maPalette}}. defaults to "black"}
  \item{byclass}{logical indicating whether samples should be clustered within each level of the phenotype or overall. defaults to F}
  \item{showaber}{logical indicating whether high level amplifications and homozygous deletions should be indicated on the plot. defaults to F}
  \item{amplif}{positive value that all observations equal or exceeding it  are marked by yellow dots indicating high-level changes. defaults to 1}
  \item{homdel}{negative value that all observations equal or below it  are marked by light blue dots indicating homozygous deletions. defaults to -0.75}
  \item{samplenames}{sample names}
  \item{vecchrom}{vector of chromosomal indeces to use for clustering and to display. defaults to 1:23}
  \item{titles}{plot title. defaults to "Image Plots"}
  \item{methodS}{clustering method to cluster samples. defaults to "ward"}
 % \item{dendPlot}{logical indicating whether dendogram needs to be
 %   drawn. defaults to T.}
  \item{categoricalPheno}{logical indicating whether phenotype is
    categorical. Continuous phenotypes are treated as "no groups" except
    that their values are dispalyed.defaults to TRUE.}
    \item{CENTROMERE}{logical indicating whether to plot the centromere
  }
}
\details{
  This functions is a more flexible version of the
  \code{\link[stats:heatmap]{heatmap}}. It can cluster within levels of categorical
  phenotype as well as all of the samples while displaying phenotype
  levels in different colors. It also uses any combination of
  chromosomes that is requested and clusters samples based on these
  chromosomes only. It draws the chromosomal boundaries and displays
  high level changes and homozygous deletions. If phenotype if not
  categorical, its values may still be displayed but groups are not formed
  and byclass = F.
  Image plot has the samples reordered according to clustering order.
}
% \value{
% }
%\references{ ~put references to the literature/web site here ~ }
\seealso{
  \code{\link[stats:heatmap]{heatmap}}
}

\keyword{hplot}
\keyword{cluster}
