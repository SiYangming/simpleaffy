\name{hmap.eset}
\alias{hmap.eset}
\title{ Draw a heatmap from an AffyBatch object}
\description{
   Given either an \code{AffyBatch} draw a heatmap. 
}
\usage{
hmap.eset(x,probesets,samples=1:length(sampleNames(x)),scluster=standard.pearson,pcluster=standard.pearson,slabs=sampleNames(x)[samples],plabs,col="bwr",min.val=NULL ,max.val=NULL,scale=FALSE,spread=6,by.fc=F,sdev=NULL,show.legend=T,title=NULL,cex=0.5)
}

\arguments{
  \item{x}{ The \code{AffyBatch} object to get the expression data from }
  \item{probesets}{ What probesets to plot, defaults to all of them }
  \item{samples}{ Which samples to plot }
  \item{scluster}{ The function to use to cluster the samples by. Can also be a dendrogram object. }
  \item{pcluster}{ The function to use to cluster the probesets by. Can also be a dendrogram object. }
  \item{slabs}{ Labels for the sample axis }
  \item{plabs}{ Labels for the probeset axis defaults to geneNames(x) }
  \item{col}{ Vector of colour values to use (see below) }
  \item{min.val}{ The minimum intensity to plot }
  \item{max.val}{ The maximum intensity to plot }
  \item{scale}{Scale each gene's clouring based on standard deviation (See below)}
  \item{spread}{If the data is scaled, how many standard deviations (or fold changes) either way should we show. If no scaling, then does nothing }
   \item{by.fc}{If the data is scaled, scale by s.d. or by fold.change? }
  \item{sdev}{A vector of standard deviaitions for each gene to be plotted. If no value is supplied these are worked out from the data.}
  \item{show.legend}{Draw a scale on the graph and show the title if supplied}
  \item{title}{The title of the graph}
  \item{cex}{Character expansion}
}

\details{
Takes an \code{AffyBatch} object and plots a heatmap. At its simplest, all that is required is an AffyBatch object (as calculated by \code{\link{call.exprs}}) and a vector supplying the probesets to plot. These can be specified by name, as an integer index or as a vector of TRUEs and FALSES. The function will try to do something sensible with the labels. If it fails you will need to specify this with plabs. The function will then draw a heatmap, coloured blue-white-red in increasing intensity, scaled so that 100% blue and 100% red appear somewhere on the image. The legend shows how these colours translate into intensity. 

Col can be used to change the colouring. "bwr" specifies blue-white-red, "rbg" specifies red-black-green, and "ryw" specifies red-yellow-white. Alternatively, a vector of arbitrary colours can be supplied (try \code{rainbow(21)}, for example).

The clustering method can also be changed by supplying, either, a function that takes a matrix of expression values and returns an \code{hclust} or \code{dendrogram} object, or alternatively, an \code{hclust} or \code{dendrogram} object itself. Setting either of these to NULL will stop the heatmap being clustered on that axis.

Scaling is somewhat more complex. If scale is TRUE, then each gene is coloured independently, on a scale based on its standard deviation. By default this is calculated for the samples that are being plotted, unless a value is supplied for sdev -- in which case this should be a vector of standard deviations, one for each probeset being plotted (and in the same order). This scaling is done after the clustering. For more details on how all of this works see the website \url{http://bioinf.picr.man.ac.uk/simpleaffy} and also look at \code{\link{hmap.pc}} which uses the scaling to plot transcripts identified as being differentially expressed.
 }

\value{ Returns an (invisible) list containing the dendrograms used for samples and probesets}

\author{ Crispin J Miller }

\seealso{ \code{\link{hmap.pc}} \code{\link{blue.white.red.cols}} \code{\link{standard.pearson}}}
\examples{
  \dontrun{ 
    eset.mas <- call.exprs(eset,"mas5")
    hmap.eset(eset.mas,1:100,1:6,col="rbg")
  }
}
\keyword{ misc }
