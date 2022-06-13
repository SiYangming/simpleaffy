\name{plot.pairwise.comparison}
\alias{plot.pairwise.comparison}
\alias{plot,PairComp}
\alias{plot,PairComp-method}

\title{ Plots a PairComp object }
\description{
   Draws a scatter plot between means from a pairwise comparison. Colours according to PMA calls and identifies 'signficant' genes yielded by a filtering
}
\usage{
\method{plot}{pairwise.comparison}(x,y=NULL,labels=colnames(means(x)),showPMA=TRUE,type="scatter",...)
}

\arguments{
  \item{x}{ A \code{PairComp} object }
  \item{y}{ A \code{PairComp} object }
  \item{labels}{ A list containing x and y axis labels }
  \item{showPMA}{ True if PMA calls are to be identified }
  \item{type}{ Can be 'scatter', 'ma' or 'volcano' }
  \item{...}{ Additional arguments to plot }
}
\details{

Takes a PairComp object (as produced by \code{pairwise.comparison}
and plots a scatter plot between the sample means. If PMA calls are
present in the \code{calls} slot of the object then it uses them to
colour the points. Present on all arrays: red; absent on all arrays:
yellow; present in all some arrays; orange. In
addition, if a second \code{PairComp} object is supplied, it
identifies spots in that object, by drawing them as black
circles. This allows, for example, the results of a
\code{pairwise.filter} to be plotted on the same graph.

If type is 'scatter' does a simple scatter plot.
If type is 'volcano' does a volcano plot.
If type is 'ma' does an MA plot.
}
\author{ Crispin J Miller }

\seealso{ \code{\link{pairwise.comparison}} \code{\link{pairwise.filter}} \code{\link{trad.scatter.plot}}}
\examples{
  \dontrun{
    pc <- pairwise.comparison(eset.mas,group="group",members=c("a","b"),spots=eset)
    pf <- pairwise.filter(pc)
    plot(pc,pf)
  }
}


\keyword{ misc }
