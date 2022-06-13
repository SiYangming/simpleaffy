\name{qc}
\alias{qc}
\alias{qc,AffyBatch-method}

\title{ Generate QC stats from an AffyBatch object }

\description{
Generate QC metrix for Affymetrix data.
}

\details{
Affymetrix recommend a series of QC metrics that should be used to
check that arrays have hybridised correctly and that sample quality is
acceptable. These are discussed in the document 'QC and Affymetrix
data' accompanying this package, and on the web at
http://bioinformatics.picr.man.ac.uk. They are described in detail in
the 'Expression Analysis Fundamentals' manual available from
Affymetrix.

Before using this function you are strongly encouraged to read the 'QC
and Affymetrix data' document, which contains detailed examples.

This function takes an \code{\link[affy:AffyBatch-class]{AffyBatch}}
object and generates a \code{\link[simpleaffy:QCStats-class]{QCStats}} 
object containing a set of QC metrics. See
\code{\link[simpleaffy]{qc.affy}} for more details. 
}
\usage{
qc(unnormalised, ...)
}
\arguments{
  \item{unnormalised}{ An AffyBatch object with nowt done to it }
  \item{\dots}{ Any other parameters }
}

\author{ Crispin J Miller }

\seealso{
  \code{\link[simpleaffy]{qc.affy}}
  \code{\link[simpleaffy]{setQCEnvironment}}
}
\examples{
  \dontrun{
    qcs <- qc(eset,eset.mas)
  }
  data(qcs)
  ratios(qcs)
  avbg(qcs)
  maxbg(qcs)
  minbg(qcs)
  spikeInProbes(qcs)
  qcProbes(qcs)
  percent.present(qcs)
  plot(qcs)
  sfs(qcs)
  target(qcs)
  ratios(qcs)
}


\keyword{ misc }
