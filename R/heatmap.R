library("methods")


standard.pearson <- function(x) {
  hclust(as.dist(1-cor(x,method="pearson")))
}

blue.white.red.cols  <- c(hsv(2/3,(10:0)/10,1),hsv(0,(1:10)/10,1))
red.black.green.cols <- c(hsv(0,1,(10:0)/10),hsv(1/3,1,(1:10)/10))
red.yellow.white.cols <- heat.colors(21)


hmap.eset <- function(x,probesets,samples=1:length(sampleNames(x)),scluster=standard.pearson,pcluster=standard.pearson,slabs=sampleNames(x)[samples],plabs,col="bwr",min.val=NULL ,max.val=NULL,scale=FALSE,spread=6,by.fc=F,sdev=NULL,show.legend=T,title=NULL,cex=0.5) {
  # set up the graphics
  par(mai=c(0,0,0,0))
  if(show.legend) {
    m <- matrix(c(0,0,1,0,0,0,0,6,2,3,4,7,8,9,0,0,5,0,0,0,0),ncol=3,nrow=7)
    layout(m,  c(1, 5, 1),c(0.33,1, 5, 1,0.1,0.2,0.1) )
  }
  else {
    m <- matrix(c(0,1,0,2,3,4,0,5,0),ncol=3,nrow=3)
    layout(m,  c(1, 5, 0.75),c(1, 5, 0.75) )
  }
  if(length(col) == 1) {
    if(col=="bwr") {
      col<-blue.white.red.cols
    }
    else {  
    if(col=="rbg") {
      col<-red.black.green.cols
    } else {
    if(col=="ryw") {
      col<-red.yellow.white.cols
    }
    }}
  }

  if(missing(probesets)) {
    probesets <- seq_along(featureNames(x))
    probesets.was.set <- T
  }
  else {
    probesets.was.set <- F
  }
  if(missing(plabs)) {
    if(probesets.was.set | is.logical(probesets) | is.integer(probesets)) {
      plabs <- featureNames(x)[probesets]
    }
    else {
      plabs <- probesets
    }    
  }

  todo <- exprs(x)[probesets,samples]

  if(!is.null(pcluster)) {
    if(is.integer(pcluster)) {
	po <- pcluster
        plot(0,0,type="n",xaxs="i",xaxt="n",yaxt="n",bty="n")
    }
    else {
      if(class(pcluster) == "dendrogram") { 
        pd <- pcluster 
      }
      else {
        pd<-as.dendrogram(pcluster(t(todo)))
      }
      po <- order.dendrogram(pd)
      plot(pd,horiz=TRUE,leaflab="none",axes=FALSE,yaxs="i")
    }
  }
  else {
    po <- 1:length(plabs)
    plot(0,0,type="n",xaxs="i",xaxt="n",yaxt="n",bty="n")
    pd <- NULL
  }

  if(!is.null(scluster)) {
    if(is.integer(scluster)) {
	so <- scluster
        plot(0,0,type="n",xaxs="i",xaxt="n",yaxt="n",bty="n")
    }
    else {
      if(class(scluster) == "dendrogram") { 
        sd <- scluster 
      }
      else {
        sd<-as.dendrogram(scluster(todo))
      }
      so <- order.dendrogram(sd)
      plot(sd,leaflab="none",axes=FALSE,xaxs="i")
    }
  }
  else {
    so <- 1:length(slabs)
    plot(0,0,type="n",xaxs="i",xaxt="n",yaxt="n",bty="n")
    sd <- NULL
  }
  


  ls <- length(so)
  lp <- length(po)

  if(!is.null("min.val")) { min.val<-min(todo) }
  if(!is.null("max.val")) { max.val<-max(todo) }

  legend.vals <- 1:length(col)* (max.val-min.val)/length(col) + min.val;

  #scale the colors if requested


  if(scale & !by.fc) {
    i <- 1;
    todo <- t(apply(2^todo,1,function(x) { 
                        minv <- min(x)
                        if(is.null(sdev)) {
                          this.sd <- sqrt(var(x))
                        }
                        else {
                          this.sd <- sdev[i];
                          i<-i + 1;
                        }
                        x <- x-minv
                        x <- x/this.sd
                        x <- x - mean(x)
                        return(x/spread)
                        }))
    
     min.val<- -spread 
     max.val<- spread
     todo[todo<min.val] <- min.val
     todo[todo>max.val] <- max.val

  }

  if(scale & by.fc) {
    todo <- t(apply(todo,1,function(x) { 
                        return(log2(2^x/mean(2^x)))
                        }))
    
     min.val<- -spread 
     max.val<- spread
     todo[todo<min.val] <- min.val
     todo[todo>max.val] <- max.val

  }

  # reorder the expression array to fit the dendogram

  todo <- todo[po,]
  todo <- todo[,so]
 
  xp <- (1:ls+1)-0.5
  yp <- (1:lp+1)-0.5

  # HURRAH! The heatmap itself!



  image(xp,yp,z=t(todo),xaxt="n",yaxt="n",bty="o",col=col,zlim=range(min.val,max.val))

  # sample labels

  plot(0,0,xlim=range(0,ls),ylim=range(0,1),type="n",xaxs="i",xaxt="n",yaxt="n",bty="n")
  text((1:ls)-0.5,1,labels=slabs[so],srt=90,adj=1,cex=cex)



  # probeset labels

  plot(0,0,xlim=range(0,1),ylim=range(0,lp),type="n",yaxs="i",xaxt="n",yaxt="n",bty="n")
  text(0,(1:lp)-0.5,labels=plabs[po],adj=0,cex=cex)


  if(is.null(title)) { title="" }

  if(show.legend) {
    plot(0,0,xlim=range(0,1),ylim=range(0,1),type="n",yaxs="i",xaxt="n",yaxt="n",bty="n")
    text(0.5,0.5,labels=title,adj=1,cex=2)

    plot(0,0,type="n",xaxt="n",xaxs="i",yaxt="n",yaxs="i",xlim=c(1,length(col)),ylim=c(0,1),bty="n")
  if(scale) {
    text(1,0.5,label=round(-spread,2),pos=4)
    text(ceiling(length(col)/4),0.5,round(-spread/2,2))
    text(ceiling(3*length(col)/4),0.5,round(spread/2,2))
    text(length(col),0.5,label=round(spread,2),pos=2)
    if(!by.fc) {
      text(ceiling(length(col)/2),0.5,label="s.d. from mean")
    }
    else {
      text(ceiling(length(col)/2),0.5,label="fold change from mean")
    }
  }
  else {
    text(1,0.5,round(legend.vals[1],2),pos=4)
    text(ceiling(length(col)/4),0.5,round(legend.vals[ceiling(length(col))/4],2))
    text(ceiling(3*length(col)/4),0.5,round(legend.vals[ceiling(3*length(col))/4],2))
    text(length(col),0.5,round(legend.vals[length(col)],2),pos=2)
    text(ceiling(length(col)/2),0.5,round(legend.vals[ceiling(length(col))/2],2))
  }

    image(x=1:length(col),y=1,z=matrix(legend.vals,nrow=length(col)),col=col,yaxt="n",ylab="",xlab="",xaxt="n",bty="n")

  }
  invisible(list(probesets=pd,samples=sd))
}

hmap.pc <- function(x,eset,samples=rownames(pData(x)),scluster=standard.pearson,pcluster=standard.pearson,slabs,plabs,col="rbg",scale=T,spread=10,by.fc=F,gp=group(x),mbrs=members(x),show.legend=T,title=NULL,cex=.10) {
  pns <- names(fc(x))
   if(length(samples)==1) {
     if(samples == "all") { sns <- rownames(pData(eset)) }
     else { sns <- samples }
   }
   else { sns <- samples }
  if(length(sns) == 0) {stop("Error 'x' contains 0 probesets")}
  if(length(sns) == 0) {stop("Error 'x' contains 0 samples")}
  if(missing(slabs)) { slabs=sns} 
  if(missing(plabs)) { plabs=pns} 
  if(scale & !by.fc) {

   nms <- names(fc(x))
   sm <- rep(0,length(nms))
   if(length(mbrs)==1) {
     if(mbrs == "all") { mbrs <- unique(as.character(pData(eset)$gp)) }
   }
   for(i in mbrs) {
     print(paste("Getting variance for member '",i,"' in group '",gp,"'.",sep=""))
     sset <- get.array.subset(eset,gp,i)
     sm <- sm + apply(exprs(sset)[nms,],1,function(y) { sqrt(var(y)) })
   }
   sm <- sm / length(mbrs)
  }
  hmap.eset(eset,samples=sns,probesets=pns,scluster=scluster,pcluster=pcluster,col=col,slabs=slabs,plabs=plabs,scale=scale,spread=spread,by.fc=by.fc,sdev=sm,show.legend=T,title=title,cex=cex)
}

