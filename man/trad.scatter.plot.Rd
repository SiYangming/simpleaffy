\name{trad.scatter.plot}
\alias{trad.scatter.plot}

\title{ Does a Traditional Scatter Plot of Expression Data }
\description{
  Plots expression data as a scatter plot with optional fold-change lines
}
\usage{
trad.scatter.plot(x, y, add = FALSE, fc.lines = log2(c(2, 4, 6, 8)), draw.fc.lines = TRUE, draw.fc.line.labels = TRUE, fc.line.col = "lightgrey", pch = 20,xlim=NULL,ylim=NULL, ...)
}

\arguments{
  \item{x}{ x coords }
  \item{y}{ y coords }
  \item{add}{ add this data to an existing graph }
  \item{fc.lines}{ Vector of intervals at which to draw fold-change lines }
  \item{draw.fc.lines}{ Draw fold change lines? }
  \item{draw.fc.line.labels}{ Label the fold change lines with the fold changes they represent? }
  \item{fc.line.col}{ The colour to draw fold change lines }
  \item{pch}{ Plotting character to use for the scatter data (see \code{plot} for more details) }
  \item{xlim}{ Range for the xaxis }  
  \item{ylim}{ Range for the yaxis }
  \item{\dots}{ Additional parameters to pass through to the underlying \code{plot} function }
}
\references{ http://bioinformatics.picr.man.ac.uk/ }
\author{ Crispin J Miller }

\seealso{ \code{\link{plot}}}


\examples{
  \dontrun{
    trad.scatter.plot(exprs(eset.rma)[,1],exprs(eset.rma)[,4])
  }
}

\keyword{ misc }
