\name{qc.affy}
\alias{qc.affy}

\title{ Generate Affymetrix Style QC Statistics }
\description{
Generate QC data for Affymetrix arrays
}


\details{
Affymetrix recommend a series of QC metrics that should be used to
check that arrays have hybridised correctly and that sample quality is
acceptable. These are discussed in the document 'QC and Affymetrix
data' accompanying this package, and on the web at
http://bioinformatics.picr.man.ac.uk. They are described in detail in
the 'Expression Analysis Fundamentals' manual available from
Affymetrix.

This function takes an (unnormalised) AffyBatch object, and (optionally) 
an ExprSet, with MAS expression calls produced by \code{\link{call.exprs}} and
generates QC metrics. If the MAS calls are not supplied these are claculated internally.

}
\usage{
qc.affy(unnormalised,normalised=NULL,tau=0.015,logged=TRUE,
	cdfn=cdfName(unnormalised))
}

\arguments{
  \item{unnormalised}{ An unnormalised raw \code{AffyBatch} object to call qc stats on }
  \item{normalised}{ The same one, processed using \code{\link{justMAS}} (contains scale factors etc.). If not supplied, then the object gets calculated internally. }
  \item{tau}{used by detection p value}
  \item{logged}{ True if the data is logged }
  \item{cdfn}{The cdf name for the array - can be used to specify a different set of probes to the default}  
}

\value{
  A QCStats object describing the supplied
  \code{\link[affy:AffyBatch-class]{AffyBatch}} 
}

\author{ Crispin J Miller }

\examples{
  \dontrun{
     qcs <- qc(eset)
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
