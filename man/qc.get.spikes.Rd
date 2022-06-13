\name{  qc.get.spikes}
\alias{  qc.get.spikes}
\alias{  qc.get.spike}
\alias{  qc.add.spike}


\title{ Retrieve QC spike probeset names for the current array type }

\description{
   Get the names of spike probesets for bioB, bioC, etc. ratios for the
   current array type. \code{\link[simpleaffy]{qc.get.probes}} is used to define the
   3'/5' ratio probesets
 }
\usage{
 qc.get.spikes()
 qc.get.spike(name)
 qc.add.spike(name,probeset)
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
  \code{\link[simpleaffy]{setQCEnvironment}}   \code{\link[simpleaffy]{qc.get.probes}}
}
\examples{
  qc.get.spikes()
  qc.add.spike("my.name","a.probesetid_at")
  qc.add.spike("another.name","another.probesetid_at")
  qc.get.spikes()
}
\keyword{ misc }
