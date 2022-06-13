\name{ qc.have.params }
\alias{ qc.have.params }

\title{ Does simpleaffy have a QC definition file for the specified array? }
\description{
  Simpleaffy requires a definition file describing the qc probes,
  spikes, alpha values, etc. for the array of interest. This is used to
  initialize the QC environment for the array (usually implicitly within
  the \code{\link[simpleaffy]{qc}} function), by a call to  
  \code{\link[simpleaffy]{setQCEnvironment}}. This function can be used
  to see if the specified array has a definition file.
}
\usage{
  qc.have.params(name)
}

\arguments{
 \item{name}{The 'clean' CDF name of the array (i.e. the result of calling
   \code{\link[affy]{cleancdfname}} on the \code{cdfName} of the
 \code{AffyBatch} object containing the  array data of interest.}
}

\value{
  True or False
}

\references{ http://bioinformatics.picr.man.ac.uk/ }
\author{ Crispin J Miller }

\seealso{
  \code{\link[simpleaffy]{setQCEnvironment}}, \code{\link[simpleaffy]{qc}},
  \code{\link[simpleaffy]{qc.ok}}, \code{\link[affy:AffyBatch-class]{cdfName}},
  \code{\link[affy]{cleancdfname}}
}
\examples{
  qc.have.params("nosucharraycdf")
  qc.have.params("hgu133plus2cdf")
  setQCEnvironment("hgu133plus2cdf")
  qc.have.params(cleancdfname("HG-U133_Plus_2"))
}
\keyword{ misc }
