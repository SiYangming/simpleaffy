\name{qc.get.array}
\alias{qc.get.array}
\alias{qc.set.array}

\title{ Get or set the name of the array for which the current QC
  environment is valid. Currently not used for anything important; is a free
  text identifier.} 

\description{
  The array name is simply a free text name for the array of interest.
}
\usage{
 qc.get.array()
 qc.set.array(name)
}

\arguments{
  \item{name}{ a free text name for the array of interest }
}

\value{
  a string
}

\references{\url{http://bioinformatics.picr.man.ac.uk/}}
\author{ Crispin J Miller }

\seealso{
  \code{\link[simpleaffy]{setQCEnvironment}}
}
\examples{
  qc.set.array("plus2")
  qc.get.array()
}
\keyword{ misc }
