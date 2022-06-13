\name{ qc.ok }
\alias{ qc.ok }

\title{ Has simpleaffy's QC environment been set up? }
\description{
  Simpleaffy requires a definition file describing the qc probes,
  spikes, alpha values, etc. for the array of interest. This is used to
  initialize the QC environment for the array (usually implicitly within
  the \code{\link[simpleaffy]{qc}} function), by a call to  
  \code{\link[simpleaffy]{setQCEnvironment}}. This function can be used
  to check if the qc environment has been set up for the current session
}
\usage{
  qc.ok()
}


\value{
  True or False
}

\references{ http://bioinformatics.picr.man.ac.uk/ }
\author{ Crispin J Miller }

\seealso{
  \code{\link[simpleaffy]{setQCEnvironment}} \code{\link[simpleaffy]{qc}}
  \code{\link[simpleaffy]{qc.have.params}}
  \code{\link[affy:AffyBatch-class]{cdfName}} 
}
\examples{
  qc.ok()
  setQCEnvironment("hgu133plus2cdf")
  qc.ok()
}
\keyword{ misc }
