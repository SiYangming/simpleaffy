\name{call.exprs}
\alias{call.exprs}
\title{ Generate Expression Summaries for Affymetrix Data }
\description{
  Generates expression summaries and normalizes Affymetrix data using either MAS5.0, GCRMA or RMA algorithms.
}
\usage{
call.exprs(x, algorithm = "rma", do.log = TRUE, sc = 100, method = NA)
}

\arguments{
  \item{x}{ an \code{AffyBatch} object }
  \item{algorithm}{ one of '"rma","rma-R","gcrma", "mas5",
  "mas5-R". "rma" and "mas5" make use of a native C-library and are
  faster than "rma-R" and "mas5-R". } 
  \item{do.log}{ return logged data if true } 
  \item{sc}{ if the mas5 algorithm is being used, sets the target
  intensity to which the chips should be scaled. } 
  \item{method}{ The algorithm used to normalise the data. Has no effect
  for "rma", defaults to quantile normalisation for "rma" and no
  normalisation for "mas5" } 
}
\value{
  An AffyBatch object containing expression summaries.
}
\references{\url{http://bioinformatics.picr.man.ac.uk/}}
\author{ Crispin J Miller }

\seealso{ \code{\link{read.affy}}, \code{\link[affy]{expresso}},
  \code{\link[affy:justrma]{justRMA}}, \code{\link{justMAS}} } 
\keyword{ misc }
\examples{ 
  \dontrun{
    eset.rma  <- call.exprs(eset,"rma"); 
    eset.mas5 <- call.exprs(eset,"mas5"); 
  }
}
