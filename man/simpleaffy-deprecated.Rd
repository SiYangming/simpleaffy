\name{ simpleaffy-deprecated }
\alias{ simpleaffy-deprecated }
\alias{getTao}
\alias{getAlpha1}
\alias{getAlpha2}
\alias{getActin3}
\alias{getActinM}
\alias{getActin5}
\alias{getGapdh3}
\alias{getGapdhM}
\alias{getGapdh5}
\alias{getAllQCProbes}
\alias{getBioB}
\alias{getBioC}
\alias{getBioD}
\alias{getCreX}
\alias{getAllSpikeProbes}
\alias{haveQCParams}

\title{ Does simpleaffy have a QC definition file for the specified array? }
\description{
  The underlying implementation of simpleaffy has changed significantly
  and it now represents QC parameters differently. In particular, it
  loads only the QC data for the specified array type. A call to any of these
  functions loads the appropriate environment specifed by \code{name}.
  They therefore been deprecated and WILL disappear
  from simpleaffy in the future.
}
\usage{
  getTao(name)
  getAlpha1(name)
  getAlpha2(name)
  getActin3(name)
  getActinM(name)
  getActin5(name)
  getGapdh3(name)
  getGapdhM(name)
  getGapdh5(name)
  getAllQCProbes(name)
  getBioB(name)
  getBioC(name)
  getBioD(name)
  getCreX(name)
  getAllSpikeProbes(name)
  haveQCParams(name)
}

\arguments{
  \item{name}{The 'clean' CDF name of the array (i.e. the result of calling
   \code{\link[affy]{cleancdfname}} on the
  \code{\link[affy:AffyBatch-class]{cdfName}} 
   of the AffyBatch object containing the
   array data of interest.}
}

\details{
  Each of these functions has been replaced by a new function of the
  form \code{qc.get.}. In order to support ratios other than gapdh and
  beta-actin, the appropriate way to get ratios is now to use
  \code{\link[simpleaffy]{qc.get.ratios}}, which will return a table
  containing all suggested ratio calculations for the
  array.
  Similarly, \code{\link[simpleaffy]{qc.get.spikes}} will return a table
  containing all spike probesets for the array.
}

\value{
  None.
}

\references{ http://bioinformatics.picr.man.ac.uk/ }
\author{ Crispin J Miller }

\seealso{
  \code{\link[simpleaffy]{setQCEnvironment}}
  \code{\link[simpleaffy]{qc}}
  \code{\link[simpleaffy]{qc.ok}}
  \code{\link[affy:AffyBatch-class]{cdfName}}
  \code{\link[affy]{cleancdfname}}
  \code{\link[simpleaffy]{qc.get.ratios}}
  \code{\link[simpleaffy]{qc.get.spikes}}
  \code{\link[simpleaffy]{qc.get.probes}}
}
\examples{
  #old
  getBioB("hgu133plus2cdf")
  getActin3("hgu133plus2cdf")
  getActinM("hgu133plus2cdf")
  getActin5("hgu133plus2cdf")
  #new
  setQCEnvironment("hgu133plus2cdf")
  qc.get.spikes()["bioB"]
  r <- qc.get.probes()
  r["actin3"]
  r["actinM"]
  r["actin5"]
}
\keyword{ misc }
