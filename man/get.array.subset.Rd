\name{get.array.subset}
\alias{get.array.subset}
\alias{get.array.subset,AffyBatch-method}
\alias{get.array.subset,ExpressionSet-method}

\title{Get a subset of arrays from an affybatch object, split by phnotypic data}
\description{
  Looks at a factor in the phenotypic data for an \code{AffyBatch} or \code{ExpressionSet} object and
  uses it to select a subset of arrays, as defined by 'members'.
}


\usage{
get.array.subset(x,group,members)
}

\arguments{
  \item{x}{An \code{ExpressionSet} or \code{AffyBatch} object.}
  \item{group}{The name of the pData column to use.}
  \item{members}{The labels within the pData column to match against.}
}

\author{Crispin J Miller}

\seealso{\code{\link{get.array.subset.affybatch}}  \code{\link{get.array.subset.exprset}}}

\examples{
  \dontrun{
    subset1 <- get.array.subset.affybatch(eset.rma,"group","A")
    subset2 <- get.array.subset.exprset(eset.rma,"group",c("A","P"))
    subset3 <- get.array.subset(eset.rma,"group","A")
  }
}

\keyword{misc}
