\name{hmap.pc}
\alias{hmap.pc}
\title{ Draw a heatmap from an PairComp object}
\description{
   Given either a \code{PairComp} object draw a heatmap. 
}
\usage{
hmap.pc(x,eset,samples=rownames(pData(x)),scluster=standard.pearson,pcluster=standard.pearson,slabs,plabs,col="rbg",scale=T,spread=10,by.fc=F,gp=group(x),mbrs=members(x),show.legend=T,title=NULL,cex=0.1)
}

\arguments{
  \item{x}{ The \code{PairComp} object to get the probeset list (and other data) from}
  \item{eset}{ The \code{AffyBatch} object containing expression data }
  \item{samples}{ Which samples to plot -- defaults to those used to calculate 'x', but can be any of the samples in eset }
  \item{scluster}{ The function to use to cluster the samples by. Can also be a dendrogram object. }
  \item{pcluster}{ The function to use to cluster the probesets by. Can also be a dendrogram object. }
  \item{slabs}{ Labels for the sample axis }
  \item{plabs}{ Labels for the probeset axis }
  \item{col}{ Vector of colour values to use (see below) }
  \item{scale}{Scale each gene's clouring based on standard deviation (See below)}
  \item{spread}{If the data is scaled, how many standard deviations (or fold changes) either way should we show. If no scaling, then does nothing }
  \item{by.fc}{If the data is scaled, do it by fold change? }
  \item{gp}{ The column in the expression set's pData object used to select the samples to plot. By default this is the one used to calculate x. }
  \item{mbrs}{The members of the 'group' column that we wish to plot. By default these are the pair used to calculate x. If 'all' is supplied then all samples are used.}
  \item{show.legend}{Draw a scale on the graph and show the title if supplied}
  \item{title}{The title of the graph}
  \item{cex}{Character expansion}
}

\details{
Takes a \code{PairComp} object and an \code{AffyBatch} object and plots a heatmap. At its simplest, all that is required are these two objects. The function will then draw a heatmap, coloured red-black-green in increasing intensity, scaled for each gene based on standard deviation. The legend shows how these colours translate into intensity. 

Col can be used to change the colouring. "bwr" specifies blue-white-red, "rbg" specifies red-black-green, and "ryw" specifies red-yellow-white. Alternatively, a vector of arbitrary colours can be supplied (try \code{rainbow(21)}, for example).

Scaling is somewhat complex. If scale is TRUE, then each gene is coloured independently, on a scale based on its standard deviation. This is calculated as follows: 'group' supplies a column in the pData object of 'eset' that is used to collect samples together (generally as replicate groups). 'members' supplies the entries within this column that are to be used. (Unless specified, the function uses the same value for 'group' and 'members' used to calculate the PairComp object). The function uses these data to calculate the standard deviation for each probeset within each set of replicates, and then calculates the average sd for each gene. This is then used to scale the data so that each probeset is plotted on a scale that shows the number of standard deviations away from the mean it is for that sample. For more details on how all of this works see the website \url{http://bioinf.picr.man.ac.uk/simpleaffy}.

Alternatively, by setting by.fc to FALSE, scaling can be done simply in terms of fold-change, in which case, spread defines the maximum and minimum fold changes to show.
 }

\value{ Returns an (invisible) list containing the dendrograms used for samples and probesets}


\author{ Crispin J Miller }

\seealso{ \code{\link{hmap.eset}} \code{\link{blue.white.red.cols}} \code{\link{standard.pearson}}}
\examples{
  \dontrun{
    pc <- pairwise.comparison(eset.mas,group="group",members=c("a","b"),spots=eset)
    pf <- pairwise.filter(pc)
    hmap.pc(pf,eset.mas)
  }
}
\keyword{ misc }
