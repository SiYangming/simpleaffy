\name{detection.p.val}
\alias{detection.p.val}

\title{ Calculate Detection p-values }
\description{
  Calculate MAS5 detection pvalues and Present Marginal Absent calls. This is an implementation based on the algorithm described in 
  Liu, Mei et al. (2002) 'Analysis of high density expression microarrays with signed-rank call algorithms', Bioinformatics 18(12) pp1593-1599.
}
\usage{
detection.p.val(x, tau = NULL,calls=TRUE,alpha1=NULL,alpha2=NULL,ignore.saturated=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ An unnormalised AffyBatch object }
  \item{tau}{ Errrmmm... tau }
  \item{alpha1}{ Present-Marginal threshold }
  \item{alpha2}{ Marginal-Absent threshold }
  \item{calls}{ if true, generate PMA calls }
  \item{ignore.saturated}{ if true do the saturation correction described in the paper, with a saturation level of 46000 }
}
\value{
  A list:
  \item{pval}{A matrix of detection p values}
  \item{call}{A matrix of PMA calls}
}
\references{ http://bioinformatics.picr.man.ac.uk/ }
\author{ Crispin J Miller }

\note{
  alpha1 and alpha2 are parameters that change according to the chip
  type you are using. If they are not specified, the function uses the
  current QC environment to find them, and attempts to set one up if it
  is not there. This is done with an internal call to the function
  \code{\link[simpleaffy]{setQCEnvironment}}. If it is unable to find the
  appropriate config files, this will cause an error. See
  \code{\link[simpleaffy]{setQCEnvironment}} for more details.
}
\seealso{
  \code{\link[simpleaffy]{setQCEnvironment}}
}
\examples{
   \dontrun{
     dpv <- detection.p.val(eset);
   }
}
\keyword{ misc }
