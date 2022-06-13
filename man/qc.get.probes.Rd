\name{  qc.get.probes}
\alias{  qc.get.probes}
\alias{  qc.get.probe}
\alias{  qc.add.probe}


\title{ Retrieve QC probeset names for the current array type }

\description{
   Get the names of probesets used to calculate 3'/5' ratios for the
   current array type. \code{\link[simpleaffy]{qc.get.spikes}} is used to set the
   spike probe names (i.e. bioB, bioC, etc.)
}
\usage{
 qc.get.probes()
 qc.get.probe(name)
 qc.add.probe(name,probeset)
}

\arguments{
  \item{name}{ A name for the given probeset. By default, this is the
    probeset identifier}
    \item{probeset}{ A probeset ID }
}

\value{
 A character array of probeset IDs, or the requested probeset ID, as appropriate.      
}

\references{ http://bioinformatics.picr.man.ac.uk/ }
\author{ Crispin J Miller }

\seealso{
  \code{\link[simpleaffy]{setQCEnvironment}}   \code{\link[simpleaffy]{qc.get.spikes}}
}
\examples{
  setQCEnvironment("hgu133plus2cdf")
  qc.get.probes()
  qc.add.probe("my.name","a.probesetid_at")
  qc.add.probe("another.name","another.probesetid_at")
  qc.get.probes()
}
\keyword{ misc }
