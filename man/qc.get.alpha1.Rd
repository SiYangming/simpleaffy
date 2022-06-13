\name{  qc.get.alpha1}
\alias{  qc.get.alpha1}
\alias{  qc.set.alpha1}
\alias{  qc.get.alpha2}
\alias{  qc.set.alpha2}
\alias{  qc.get.tau }

\title{ Get or set the alpha values for the current QC environment }
\description{
  Alpha1 and Alpha2 are used to define the P/M/A thresholds for
  detection calling algorithm see -
  \code{\link{detection.p.val}}. These are array dependent, these
  functions set or get their values.
  Tau is a constant parameter within the calculation and is not array
  specific.
}
\usage{
 qc.get.alpha1()
 qc.set.alpha1(value)
 qc.get.alpha2()
 qc.set.alpha2(value)
 qc.get.tau()
}

\arguments{
  \item{value}{ A double representing the alpha1 or alpha2 threshold for
    defining detection calls. See
    \code{\link{detection.p.val}} for more details.}
  
}

\value{
  \code{qc.set.alpha1} and \code{qc.set.alpha2} return
  nothing. \code{qc.get.alpha1} and \code{qc.get.alpha2} return a double.
}

\references{ http://bioinformatics.picr.man.ac.uk/ }
\author{ Crispin J Miller }

\seealso{
  \code{\link{detection.p.val}}
}
\examples{
  setQCEnvironment("hgu133plus2cdf")
  qc.get.alpha1()
  qc.get.alpha2()
  qc.set.alpha1(0.05)
  qc.get.alpha1()
  qc.set.alpha2(0.05)
  qc.get.alpha2()
}
\keyword{ misc }
