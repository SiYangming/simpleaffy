"trad.scatter.plot" <-
function(x,y,add=FALSE,fc.lines=log2(c(2,4,6,8)),draw.fc.lines=TRUE,draw.fc.line.labels=TRUE,fc.line.col="lightgrey",pch=20,xlim=NULL,ylim=NULL,...) {
   if(!add) {
     mx <- max(x,y);     
     mn <- min(x,y);
     if(is.null(xlim)) { xlim=range(mn,mx) }
     if(is.null(ylim)) { ylim=range(mn,mx) }
     plot(x,y,xlim=xlim,ylim=ylim,pch=pch,...);
     if(draw.fc.lines) {
       for( lne in fc.lines) {
         xc <- xlim 
         yc <- c(ylim[1]+lne,ylim[2]+lne);
         lines(xc,yc,col=fc.line.col); 
         if(draw.fc.line.labels) {text(mn+0.25,lne+mn,2^lne,col=fc.line.col); }
         xc <- xlim
         yc <- c(ylim[1]-lne,ylim[2]-lne);
         lines(xc,yc,col=fc.line.col); 
          if(draw.fc.line.labels) {text(lne+mn,mn+0.25,2^lne,col=fc.line.col); }
       }
     }
   }
   else {
     points(x,y,pch=pch,...); 
   }
}
