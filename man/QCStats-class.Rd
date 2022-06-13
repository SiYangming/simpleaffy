\name{QCStats-class}
\docType{class}
\alias{QCStats-class}

\alias{avbg,QCStats-method}
\alias{maxbg,QCStats-method}
\alias{minbg,QCStats-method}
\alias{percent.present,QCStats-method}
\alias{sfs,QCStats-method}
\alias{target,QCStats-method}
\alias{spikeInProbes,QCStats-method}
\alias{qcProbes,QCStats-method}
\alias{ratios,QCStats-method}
\alias{arrayType,QCStats-method}

\alias{avbg-method}
\alias{maxbg-method}
\alias{minbg-method}
\alias{percent.present-method}
\alias{sfs-method}
\alias{target-method}
\alias{spikeInProbes-method}
\alias{qcProbes-method}
\alias{ratios-method}
\alias{arrayType-method}

\alias{avbg}
\alias{maxbg}
\alias{minbg}
\alias{percent.present}
\alias{sfs}
\alias{target}
\alias{spikeInProbes}
\alias{qcProbes}
\alias{ratios}
\alias{arrayType}


\alias{plot,QCStats,ANY-method}

\title{Class "QCStats" }
\description{ Holds Quality Control data for a set of Affymetrix arrays }
\section{Objects from the Class}{
Objects can be created by calls of the form
\code{qc(\link[affy:AffyBatch-class]{AffyBatch})}. 
}
\section{Slots}{
  \describe{
    \item{\code{scale.factors}:}{Object of class \code{"numeric"} Scale factors used to scale the chips to the specified target intensity }
    \item{\code{target}:}{Object of class \code{"numeric"} The target intensity to which the chips were scaled }
    \item{\code{percent.present}:}{Object of class \code{"numeric"} Number of genes called present }
    \item{\code{average.background}:}{Object of class \code{"numeric"} The average background for the arrays }
    \item{\code{minimum.background}:}{Object of class \code{"numeric"} The minimum background for the arrays }
    \item{\code{maximum.background}:}{Object of class \code{"numeric"}
	  The maximum background for the arrays }
	\item{\code{bioBCalls}:}{Object of class \code{"character"} The
	  detection PMA (present / marginal / absent) calls of bioB spike-in
	  probes }
    \item{\code{spikes}:}{Object of class \code{"list"} spiked in probes (e.g. biob, bioc...) }
    \item{\code{qc.probes}:}{Object of class \code{"list"} qc probes
      (e.g. gapdh 3,5,M,...) }
    \item{\code{arraytype}:}{The
      \code{\link[affy:AffyBatch-class]{cdfName}} of the 
      \code{\link[affy:AffyBatch-class]{AffyBatch}} object used to
      create the object } 
  }
}
\section{Methods}{
  \describe{
    \item{avbg}{\code{signature(object = "QCStats")}:  average background}
    \item{maxbg}{\code{signature(object = "QCStats")}: maximum background }
    \item{minbg}{\code{signature(object = "QCStats")}: minimum background }
    \item{spikeInProbes}{\code{signature(object = "QCStats")}: the spike-in QC probes }
    \item{qcProbes}{\code{signature(object = "QCStats")}: the gapdh and actin QC probes }
    \item{percent.present}{\code{signature(object = "QCStats")}: no probesets called present }
    \item{plot}{\code{signature(x = "QCStats")}: Plot a QCStats object }
    \item{sfs}{\code{signature(object = "QCStats")}: scale factors }
    \item{target}{\code{signature(object = "QCStats")}: target scaling }
    \item{ratios}{\code{signature(object = "QCStats")}: 5'3' and 5'M
      ratios for QC Probes }
    \item{arrayType}{\code{signature(object = "QCStats")}: The type of
      array for which this QC stats object was generated }
  }
}

\author{ Crispin J Miller }

\seealso{
  \code{\link[simpleaffy]{qc}}
}

\keyword{classes}
