\name{all.present}
\alias{all.present}

\title{ Filter by PMA call }
\description{
 must be present in at least no arrays to be called present
}
\usage{
  \method{all}{present}(x,calls,no = "all")
}

\arguments{
   \item{x}{ An object to filter }
   \item{calls}{ A matrix of PMA calls }
   \item{no}{ How many in a row to pass the filter? If 'all' then all must be present}
}

\value{
  A probesetid
}

\author{ Crispin J Miller }

\examples{
  \dontrun{
     all.present(eset,calls,dim(calls)[2])
  }
}

\keyword{ misc }
