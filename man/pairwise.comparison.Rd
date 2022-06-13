\name{pairwise.comparison}
\alias{pairwise.comparison}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Compute pairwise comparison statistics between two experimental groups }
\description{
  Generate fold changes, t-tests and means for a pair of experimental groups
}
\usage{
  pairwise.comparison(x,group,members=NULL,spots=NULL,a.order=NULL,b.order=NULL,method="unlogged",logged=TRUE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an \code{ExpressionSet} object.}
  \item{group}{column in pData(x).}
  \item{members}{labels in group.}
  \item{spots}{unnormalised AffyBatch data for this experiment - if included, results in PMA calls and detection p-values being generated}
  \item{a.order}{For a comparison with matched pairs, the ordering of the first group of replicates}
  \item{b.order}{For a comparison with matched pairs, the ordering of the second group of replicates}
  \item{method}{What method should be used to calculate the average for the fold-change - can be either "logged","unlogged","median"}
  \item{logged}{Whether the input data is logged or not}
}
\details{
  Given an ExpressionSet object, generate quick stats for pairwise
comparisons between a pair of experimental groups. If a.order and
b.order are specified then a paired sample t-test will be conducted
between the groups, with the arrays in each group sorted according to
the ordering specified. By default, the function assumes that the
expression values are logged (this can be changed with the parameter
"logged"). The fold-changes are computed from the average values
across replicates. Unless you specify otherwise, this is done using the mean of the
unlogged values (i.e. logged data is first unlogged, the mean
calculated, and the result re-logged). The parameter "method", allows
the mean of the logged values or their median to be used
instead. T-tests are always computed with the logged data.  
} 

\value{
  A Pairwise comparison object.
}

\references{http://bioinformatics.picr.man.ac.uk/}

\author{Crispin J Miller}

\examples{ 
  \dontrun{
    pc <- pairwise.comparison(eset.rma,"group",c("A","P")) 
  }
}

\keyword{misc}
