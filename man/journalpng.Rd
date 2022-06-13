\name{journalpng}
\alias{journalpng}
\alias{screenpng}

\title{ Produce a device for producing artwork for presentations and journals }
\description{
  journalpng generates a device to print a 4 x 4 inch 300 dpi figure (by default).
  screenpng does the same, but 72dpi.
}

\usage{
journalpng(file="figure.png",width=4, height=4,res=300)
screenpng(file="figure.png",width=4, height=4,res=72)
}


\arguments{
  \item{file}{ the file to write the figure to }
  \item{width}{ the width of the figure }
  \item{height}{ its height }
  \item{res}{ resolution in dots-per-inch }
}

\value{
  A table containing annotation data
}
\references{ http://bioinformatics.picr.man.ac.uk/ }
\author{ Crispin J Miller }


\examples{
   \dontrun{
     journalpng(file="results/figure1.png"); # starts a new device
     trad.scatter.plot(exprs(eset)[,1],exprs(eset)[,2])
     dev.off();    # writes the file at this point.
  }
}
\keyword{ misc }
