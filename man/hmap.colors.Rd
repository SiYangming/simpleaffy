\name{blue.white.red.cols}

\alias{blue.white.red.cols}
\alias{red.black.green.cols}
\alias{red.yellow.white.cols}


\title{ Generate colourings for heatmaps }
\description{
  Produces standard colourings for heatmaps.
}
\usage{
  blue.white.red.cols(x)
  red.black.green.cols(x)
  red.yellow.white.cols(x)
}

\arguments{
   \item{x}{ How many colours to make }
}

\value{
  A vector of colors
}

\author{ Crispin J Miller }

\examples{
  \dontrun{
       cols <- blue.white.red.cols(21)
  }
}

\seealso{
  \code{hmap}
  \code{hmap.eset}
  \code{hmap.pc}
}
\keyword{ misc }
