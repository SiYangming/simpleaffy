\name{get.fold.change.and.t.test}
\alias{get.fold.change.and.t.test}

\title{Compute fold change and t-test statistics between two experimental groups}
\description{
  Generate fold changes (and possibly means) for a pair of experimental groups
}
\usage{
  get.fold.change.and.t.test(x,group,members,logged = TRUE,a.order=NULL,b.order=NULL,method=c("unlogged","logged","median"));
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an \code{ExpressionSet} object.}
  \item{group}{column in pData(x).}
  \item{members}{labels in group.}
  \item{logged}{is the AffyBatch data logged?}
  \item{a.order}{For a pairwise comparison the ordering of the first group of replicates}
  \item{b.order}{For a pairwise comparison the ordering of the second group of replicates}
  \item{method}{What method should be used to calculate the average for the fold-change - can be either "logged","unlogged","median" }
}
\details{
  Given an ExpressionSet object, generate quick stats for pairwise
comparisons between a pair of experimental groups. If a.order and
b.order are specified then a paired sample t-test will be conducted
between the groups, with the arrays in each group sorted according to
the ordering specified.

  The fold-changes are computed from the average values across
replicates. By default this is done using the mean of the unlogged
values. The parameter, method allows the mean of the logged values or
the median to be used instead. T-tests are always computed with the
logged data.
}
\value{
  An object of class PairComp
}
\references{http://bioinformatics.picr.man.ac.uk/}
\author{Crispin J Miller}

\examples{
  \dontrun{
    pc <- get.fold.change.and.t.test(eset.rma,"group",c("A","P"))
  }
}

\keyword{misc}
