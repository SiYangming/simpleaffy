\name{all.present.in.group}
\alias{all.present.in.group}

\title{ Filter by PMA call }
\description{
   Filters an object by PMA calls. Must be called present in at leset 'no' elements in at least one of the replicate sets in the factor 'group'
}
\usage{
  \method{all}{present.in.group}(x,group,members,calls,no = "all")
}

\arguments{
   \item{x}{ An object to filter }
   \item{group}{ The factor to filter by }
   \item{members}{ The members in the group to check. If null, checks all possible ones}
   \item{calls}{ A matrix of PMA calls }
   \item{no}{ How many in a row to pass the filter? If 'all' then all must be present}
}

\value{
  A probesetid
}

\author{ Crispin J Miller }

\examples{
  \dontrun{
     all.present.in.group(eset,calls,"line",dim(calls)[2])
  }
}

\keyword{ misc }
