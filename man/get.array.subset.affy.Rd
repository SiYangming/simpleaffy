\name{get.array.subset.affybatch}
\alias{get.array.subset.affybatch}
\alias{get.array.subset.exprset}

\title{Get a subset of arrays from an affybatch object, split by phnotypic data}
\description{
  Looks at a factor in the phenotypic data for an \code{AffyBatch} or \code{ExpressionSet} object and
  uses it to select a subset of arrays, as defined by 'members'.
}
\usage{
  get.array.subset.affybatch(x, group, members)
  get.array.subset.exprset(x, group, members)
}

\arguments{
  \item{x}{An \code{AffyBatch} or \code{ExpressionSet} object.}
  \item{group}{The name of the pData column to use.}
  \item{members}{The labels within the pData column to match against.}
}

\value{
  An \code{AffyBatch} or \code{ExpressionSet} (as appropriate) containing the selected subset of chips.
}

\details{
  Subsetting an \code{AffyBatch} object by array is achieved using \code{[x,]}, while the same is achieved for an \code{ExpressionSet} by \code{[,x]}. Hence the two different functions. In general the generic method \code{\link{get.array.subset}} should be used - since it sorts this all out automatically.
}

\author{Crispin J Miller}

\examples{
  \dontrun{
    subset1 <- get.array.subset.affybatch(eset.rma,"group","A")
    subset2 <- get.array.subset.exprset(eset.rma,"group",c("A","P"))
    subset3 <- get.array.subset(eset.rma,"group","A")
  }
}
\keyword{misc}
