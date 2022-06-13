\name{get.annotation}
\alias{get.annotation}
\alias{write.annotation}
\alias{results.summary}


\title{ Get annotation data for a gene list}
\description{
  Takes a vector of probeset names and a CDF name. Produces a table of
annotations, containing gene name, description, sequence accession
number and unigene accession number for each probeset. In addition,
write.annotation is a utility function that outputs the annotation
data in a form suitable for loading into excel and results.summary
takes the outut of pairwise.comparison or pairwise.filter and spits
out a table with the means of the replicates the fold-change between
them (log2) and t-test p-values. This is followed by a table of
annotation (produced by get.annotation).
}

\usage{
get.annotation(x, cdfname,verbose=FALSE)
write.annotation(summary,file="results/annotation.table.xls")
results.summary(results,cdfname)
}

\arguments{
  \item{x}{ a vector of probe names }
  \item{cdfname}{ the name of the chip (as produced by cdfName(AffyBatch) }
  \item{verbose}{ print out information if problems are found looking things up in the annotation data }
  \item{summary}{ a table of data to write in a format appropriate to read into Excel }
  \item{file}{ a table delimited file }
  \item{results}{ a PairComp object, as produced by pairwise.comparison and pairwise.filter }
}

\value{
  A table containing annotation data
}
\references{ http://bioinformatics.picr.man.ac.uk/ }
\author{ Crispin J Miller }


\examples{
 \dontrun{
   pw       <- pairwise.comparison(eset.rma,"group",c("A","P"))
   pw.filtered <- pairwise.filter(pw)
   summary     <- results.summary(pw.filtered,"hgu133a")
   write.annotation(file="spreadsheet.xls",summary)
  }
}
\keyword{ misc }
