\name{fractionAltered}
\alias{fractionAltered}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to compute fraction of genome altered for each sample}
\description{
 This function outputs lists containing proportions of the genome that
 are gained and lost for each sample. 
}
\usage{
fractionAltered(input, thres = 0.25, factor = 2.5, 
         samplenames = colnames(input$M), 
         chrominfo = chrominfo.basepair)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{An object of type \code{\link[limma]{MAList}} or \code{\link[snapCGH:SegList]{SegList}}}
  \item{thres}{either a vector providing unique threshold for
          each sample or a vector of the same length as number of
          samples providing sample-specific
          threshold. Clone is considered to be gained if
          it is above the threshold and lost if it is below negative
          threshold. Defaults to 0.25}
  \item{factor}{specifies the number by which experimental
          variability should be multiples. Used only when tumor
          specific variability in 'input' is not NULL or when factor is greater than 0. Defaults to
          2.5.}
  \item{samplenames}{Sample names. Default is the column names of the
    \code{aCGHList} object}
%  \item{acghList}{this object should be input}
  \item{chrominfo}{A chromosomal information associated with the mapping of the
          data. Default is \code{chrominfo.basepair}}
  }
%\details{
%  ~~ If necessary, more details than the __description__  above ~~
%}
\value{
  \item{gainP}{Vector of proportion of genome gained for each sample}
  \item{lossP}{Vector of proportion of genome lost for each sample}
 
}
%\references{ ~put references to the literature/web site here ~ }
\author{Jane Fridlyand, Ritu Roydasgupta}
%\note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{\code{\link{plotSegmentationSummary}}}

\keyword{htest}

