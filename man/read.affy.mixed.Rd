\name{read.affy.mixed}
\alias{read.affy.mixed}

\title{ Read a Set of .CEL Files and Phenotypic Data from mixed chip types}
\description{
 Reads the specified file, which defines phenotypic data for a
 set of .CEL files. Reads the specified files into an
 \code{\link[affy:AffyBatch-class]{AffyBatch}} object and 
 then creates a \code{\link[Biobase]{phenoData}} object, defining the
 experimental factors  
 for those chips. This function deals with different array types by
 generating a pseudo arrayset containing only the probes in common. It
 does this by finding the smallest chip type in the set, and using this
 as a template. Probesets that aren't shared are set to 0. Other
 probesets are copied in. Note that this means that spots that were in
 one place on one array, appear to be at a different place on
 another. What this does to position specific background correction
 algorithms (such as mas5) is left as an exercise to the
 reader). Beware... 
}
\usage{
read.affy.mixed(covdesc = "covdesc",path=".", ...)
}

\arguments{
  \item{covdesc}{ 
    A white space delimited file suitable for reading as a
  \code{\link{data.frame}}. The first column (with no column name)   
     contains the names(or paths to) the .CEL files to read. Remaining
  columns (with names) represent experimental factors 
     for each chip. these become elements of the
  \code{\link[Biobase]{phenoData}} object. 
  }
  \item{\dots}{ extra functions to pass on to \code{ReadAffy} }
  \item{path}{ The path to prefix the filenames with before calling
  \code{ReadAffy} } 
}
\value{
  An AffyBatch object
}
\references{\url{ http://bioinformatics.picr.man.ac.uk/}}
\author{ Crispin J Miller }

\seealso{ \code{\link[affy:read.affybatch]{ReadAffy}},
  \code{\link[affy:AffyBatch-class]{AffyBatch}} \code{\link{data.frame}}
  \code{\link[Biobase]{phenoData}}} 
\keyword{ misc }
\examples{ 
  \dontrun{
    eset <- read.affy.mixed(); # read a set of CEL files
  
    eset.rma <- call.exprs(eset,"rma");   
  }
}
