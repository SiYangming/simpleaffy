\name{  qc.get.ratios}
\alias{  qc.get.ratios}
\alias{  qc.get.ratio}
\alias{  qc.add.ratio}


\title{ Retrieve pairs of probesets used for calculating 3'/5' ratios }

\description{
   Get the names of the qc probesets used to define the
   3'/5' ratios.
 }
\usage{
 qc.get.ratios()
 qc.get.ratio(name)
 qc.add.ratio(name,probeset1,probeset2)
}

\arguments{
  \item{name}{ A name for the given ratio calculation (such as gapdh3/5)}
  \item{probeset1}{ A probeset ID }
  \item{probeset2}{ A probeset ID }
}

\value{
 A list, each element with a name like gapdh3/5 and comprising of a
 two-element character vector of probeset IDs.      
}

\references{ http://bioinformatics.picr.man.ac.uk/ }
\author{ Crispin J Miller }

\seealso{
  \code{\link[simpleaffy]{setQCEnvironment}}   \code{\link[simpleaffy]{qc.get.probes}}
}
\examples{
  setQCEnvironment("hgu133plus2cdf")
  qc.get.ratios()
  qc.add.ratio("a.name","probeset1.id","probeset2.id")
  qc.get.ratio("a.name")
}
\keyword{ misc }
