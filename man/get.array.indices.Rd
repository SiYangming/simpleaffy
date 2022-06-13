\name{get.array.indices}
\alias{get.array.indices}
\alias{get.array.indices,AffyBatch-method}
\alias{get.array.indices,ExpressionSet-method}

\title{Find arrays in an AffyBatch object defined by their phenoData}
\description{
  Given an \code{AffyBatch} object, looks at its phenoData slot to find the factor, or column specified by 'group' and searches that column for entries supplied in 'members'. Returns the indices of these rows. For example, in a six chip AffyBatch object, x,  with a column 'treatment' containing 'c','c',t1','t2','t1','t2', a call to \code{get.array.indices(x,\"treatment\",c(\"c\",\"t1\"))} would return c(1,2,3,5).
}

\usage{
get.array.indices(x,group,members)
}

\arguments{
  \item{x}{An \code{ExpressionSet} or \code{AffyBatch} object.}
  \item{group}{The name of the pData column to use.}
  \item{members}{The labels within the pData column to match against.}
}

\author{Crispin J Miller}

\examples{
  \dontrun{
    indices3 <- get.array.indices(eset.rma,"group","A")
  }
}

\keyword{misc}
