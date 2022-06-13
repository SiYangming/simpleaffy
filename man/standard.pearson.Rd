\name{standard.pearson}
\alias{standard.pearson}

\title{ A clustering function based on pearson correlation }
\description{
  Given a matrix of values, uses hclust and cor to generate a clustering based on 1-Pearson correlation
}
\usage{
  standard.pearson(x)
}

\arguments{
   \item{x}{ A matrix of data }
}

\value{
  The result of performing an hclust
}

\author{ Crispin J Miller }

\examples{
  \dontrun{
	y <- standard.pearson(x)
  }
}

\seealso{
  \code{hmap}
  \code{hmap.eset}
  \code{hmap.pc}
}
\keyword{ misc }
