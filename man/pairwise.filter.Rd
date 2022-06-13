\name{pairwise.filter}
\alias{pairwise.filter}

\title{ Filter pairwise comparison statistics between two experimental groups }

\description{
  Given the results of a pairwise.comparison, filter the resulting gene list on expression level, PMA calls (if available), fold change and t-test statistic. 

  min.exp and min.exp.no allow the data to be filtered on intensity (where min.exp.no specifies the minimum number of arrays that must be above the threshold 'min.exp' to be allowed through the filter).
  
  PMA filtering is done when min.present.number is greater than 0.

  min.present.no allows arrays to be filtered by PMA call. A number or 'all' must be specified. If a number, then the at least this many arrays must be called present, if 'all', then all arrays must be called present. 

  present.by.group specifies whether PMA filtering is to be done on a per-group basis or for all arrays at once. If 'false' then the experiment is treated as a single group (i.e. a probeset passes the filter if it is called present on at least min.present.number arrays in the whole experiment. If 'true' then it must be called present on at least this many arrays in one or more groups. See the vignette for more details.
}

\usage{
  pairwise.filter(object,min.exp=log2(100),min.exp.no=0,min.present.no=0,present.by.group=T,fc=1.0,tt=0.001)
}


\arguments{
  \item{object}{ a 'PairComp' object }
  \item{min.exp}{ Filter genes using a minimum expression cut off }
  \item{min.exp.no}{ A gene must have an expression intensity greater than 'min.exp' in at least this number of chips }
  \item{min.present.no}{ A gene must be called present on at least this number of chips }
  \item{present.by.group}{ If true, then the probeset must be called Present on at least min.present.number arrays in any of the replicate sets used to generate the PairComp object being filtered. If false, then must be called present on at least min.present.no of the arrays in the whole experiment }
  \item{fc}{ A gene must show a log2 fold change greater than this to be called significant }
  \item{tt}{ A gene must be changing with a p-score less than this to be called significant }
}
\value{
  A 'PairComp' object filtered to contain only the genes that pass the specified filter parameters.
}
\references{ http://bioinformatics.picr.man.ac.uk/ }
\author{ Crispin J Miller }



\examples{ 
  \dontrun{
    pc <- pairwise.comparison(eset.rma,"group",c("A","P"))
    pf <- pairwise.filter(pc,tt=0.01);
  }
}

\keyword{ misc }
