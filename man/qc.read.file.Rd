\name{ qc.read.file }
\alias{ qc.read.file }

\title{ Read a file defining the QC parameters for a specified array and
set up the QC Environment }
\description{
  Affymetrix define a series of QC parameters for their arrays. Many of
  these rely on specific probeset that differ between arrays and are
  used to calculate things like 3'/5' ratios. See \code{\link[simpleaffy]{qc.affy}}
  for more details. This is usually done by a call to
  \code{\link[simpleaffy]{setQCEnvironment}}; the function described here is the one
  that does the actual loading of the
  configuration file.  See the package vignette for details of the config
  file's syntax.
}
\usage{
  qc.read.file(fn)
}

\arguments{
  \item{fn}{ full path and name of the file to load }
}

\value{
  none.
}

\references{ http://bioinformatics.picr.man.ac.uk/ }
\author{ Crispin J Miller }

\seealso{
  \code{\link[simpleaffy]{setQCEnvironment}}
}
\examples{
  fn <- system.file("extdata","hgu133plus2cdf.qcdef",package="simpleaffy")
  qc.read.file(fn)
  qc.get.spikes()
  qc.get.probes()
  qc.get.ratios()
}
\keyword{ misc }
