\name{runDNAcopy}
\alias{runDNAcopy}
\title{Results of segmenting an MAList data object using the DNAcopy library}
\description{
  The results of segmenting data from copy number array experiments from
  programs such as circular binary segmentation (CBS).  This function requires
  the library \code{DNAcopy} to be loaded.
}
\usage{
runDNAcopy(input, smooth.region=2, outlier.SD.scale = 4, smooth.SD.scale = 2, trim=0.025, alpha = 0.01, p.method = c("hybrid", 
    "perm"), kmax = 25, nmin = 200, undo.splits = c("none", "prune", "sdundo"), 
    undo.prune = 0.05, undo.SD = 3, nperm = 10000, eta = 0.05, sbdry=NULL)
}
\arguments{
  \item{input}{An object of class \code{\link[limma]{MAList}} or \code{\link[snapCGH:SegList]{SegList}}}
  \item{smooth.region}{number of points to consider on the left and the
    right of a point to detect it as an outlier.}
  \item{outlier.SD.scale}{the number of SDs away from the nearest point
    in the smoothing region to call a point an outlier.}
  \item{smooth.SD.scale}{the number of SDs from the median in the
    smoothing region where a smoothed point is positioned.}
  \item{trim}{proportion of data to be trimmed for variance calculation
    for smoothing outliers and undoing splits based on SD.}
  \item{alpha}{significance levels for the test to accept change-points.}
  \item{p.method}{method used for p-value computation.  For the "perm"
    method the p-value is based on full permutation.  For the "hybrid"
    method the maximum over the entire region is split into maximum of
    max over small segments and max over the rest.  Approximation is
    used for the larger segment max. Default is hybrid.}
  \item{kmax}{the maximum width of smaller segment for permutation
    in the hybrid method.}
  \item{nmin}{the minimum length of data for which the approximation of
    maximum statistic is used under the hybrid method.}
  \item{undo.splits}{A character string specifying how change-points are
    to be undone, if at all.  Default is "none".  Other choices are
    "prune", which uses a sum of squares criterion, and "sdundo", which 
    undoes splits that are not at least this many SDs apart.}
  \item{undo.prune}{the proportional increase in sum of squares allowed
    when eliminating splits if undo.splits="prune".}
  \item{undo.SD}{the number of SDs between means to keep a split if
    undo.splits="sdundo".} 
  \item{nperm}{number of permutations used for p-value computation.}
  \item{eta}{the probability to declare a change conditioned on the permuted statistic exceeding the observed statistic exactly j (= 1,...,nperm*alpha) times.}
  \item{sbdry}{the sequential boundary used to stop and declare a change. This boundary is a function of nperm, alpha and eta. It can be obtained using the function "getbdry" and used instead of having the "segment" function compute it every time it is called.}
}
\value{
  The function returns an object of class \code{\link[snapCGH:SegList]{SegList}}
}
\seealso{
  \code{\link[DNAcopy]{segment}}
  \code{\link[limma]{MAList}}
  \code{\link{runHomHMM}}
  \code{\link{runGLAD}}
  \code{\link[snapCGH:SegList]{SegList}}
}
\author{Mike Smith, based upon DNAcopy help files written by E. S. Venkatraman and Adam Olshen}
\keyword{methods}
