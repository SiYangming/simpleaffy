\name{bg.correct.sa}
\alias{bg.correct.sa}
\title{ Simpleaffy Implementation of Mas5 Background Correction }
\description{
  Implements the MAS5.0 background correction functions as described in Affy's
  'Statistical Algorithms Description Document'.
}
\usage{
  bg.correct.sa(unnormalised,grid=c(4,4))
}
\arguments{
  \item{unnormalised}{ An unnormalised AffyBatch object }
  \item{grid}{The dimensions of the grid to divide the chip into for background correction.}
}
\value{
  An AffyBatch object
}
\references{ \url{http://bioinformatics.picr.man.ac.uk/}
\url{http://www.affymetrix.com/support/technical/technotes/statistical_reference_guide.pdf}}
\author{ 
  Crispin J Miller 
}
\examples{ 
   \dontrun{
     eset.bg.mas <- bg.correct.sa(eset);   
   }
}
\seealso{\url{http://www.affymetrix.com/support/technical/technotes/statistical_reference_guide.pdf}}
\keyword{ misc }
