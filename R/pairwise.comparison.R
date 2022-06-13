# holds the results of a pairwise comparison
setClass("PairComp",representation(means="matrix",fc="numeric",tt="numeric",calls="matrix",group="character",members="character",pData="data.frame",calculated.from="ExpressionSet"))

#accessor methods
setGeneric("means", function(object) standardGeneric("means"))
setMethod("means","PairComp",function(object) object@means)

setGeneric("fc", function(object) standardGeneric("fc"))
setMethod("fc","PairComp",function(object) object@fc)

setGeneric("tt", function(object) standardGeneric("tt"))
setMethod("tt","PairComp",function(object) object@tt)

setGeneric("calls", function(object) standardGeneric("calls"))
setMethod("calls","PairComp",function(object) object@calls)

setGeneric("group", function(object) standardGeneric("group"))
setMethod("group","PairComp",function(object) object@group)

setGeneric("members", function(object) standardGeneric("members"))
setMethod("members","PairComp",function(object) object@members)



setMethod(Biobase::pData,"PairComp",function(object) object@pData)

setGeneric("calculated.from", function(object) standardGeneric("calculated.from"))
setMethod("calculated.from","PairComp",function(object) object@calculated.from)


##subseting. can only happen by gene

setMethod("[", "PairComp", function(x,i,j,...,drop=FALSE) {
  if(nrow(calls(x))>0) {calls <- calls(x)[i,,...,drop=FALSE];}
  else { calls <- matrix(nrow=0,ncol=0); }
  y <- new ("PairComp",means=means(x)[i,,...,drop=FALSE],fc=fc(x)[i,...,drop=FALSE],tt=tt(x)[i,...,drop=FALSE],calls=calls,group=group(x),members=members(x),pData=pData(x),calculated.from=calculated.from(x))
  return(y)
})

setReplaceMethod("[", "PairComp", function(x, i,j,...,value) {
  stop("operation not supported");
})



"get.array.subset.exprset" <-
function(x,group,members) {
  pd <- pData(x);
  grp <- pd[,colnames(pd) == group];
  return(x[,is.element(grp,members)]);
}

"get.array.subset.affybatch" <-
function(x,group,members) {
  pd <- pData(x);
  grp <- pd[,colnames(pd) == group];
  return(x[,is.element(grp,members)]);
}

"get.array.indices" <-
function(x,group,members) {
  pd <- pData(x);
  grp <- pd[,colnames(pd) == group];
  i <- 1:length(grp)
  return(i[is.element(grp,members)]);
}

setGeneric("get.array.subset", function(x,group,members) standardGeneric("get.array.subset"))
setMethod("get.array.subset","AffyBatch",get.array.subset.affybatch);
setMethod("get.array.subset","ExpressionSet",get.array.subset.exprset);
 
"get.fold.change.and.t.test" <- function(x,group,members,logged = TRUE, a.order=NULL,b.order=NULL,method=c("unlogged","logged","median")) {
  
  a.samples <- exprs(get.array.subset(x,group,members[1]));
  b.samples <- exprs(get.array.subset(x,group,members[2]));
 
  pw <- FALSE;

  if(!is.null(a.order)) { 
    a.samples <- a.samples[,a.order];
    if(!is.null(b.order)) { 
      b.samples <- b.samples[,b.order]; 
      pw <- TRUE;
 
    }
    else {
     stop("Both a.order and b.order must be specified for a paired t-test");
    }
  }
  if(length(a.order) != length(b.order)) { stop("a.order and b.order must be the same length in a paired t test") }
  method <- match.arg(method)
  m <- switch(method,
              logged   = 2,
              unlogged = 1,
              median   = 3);

  a.samples.array <- as.double(t(a.samples));
  b.samples.array <- as.double(t(b.samples));

  nacol <- as.integer(length(a.samples[1,]));
  ngene <- as.integer(length(a.samples[,1]));
  nbcol <- as.integer(length(b.samples[1,]));

  if(class(logged) != "logical") stop("Parameter 'logged' should be TRUE or FALSE")
  if((nacol == 1) | (nbcol == 1))  warning("There was only one sample in one (or both) of your sample groups. Not computing t-tests - instead, returning 0.0 for p-scores...");

  c.res <- .C("FCM",a.samples.array,b.samples.array,nacol,nbcol,ngene,as.logical(logged),as.integer(m),ma = double(ngene),mb = double(ngene),fc = double(ngene),PACKAGE="simpleaffy")
  means <- cbind(c.res$ma,c.res$mb);
  colnames(means) <- members;
  fc <- c.res$fc;

  if(!pw) {  
    ai <- dim(a.samples)[2]    
    bi <- dim(b.samples)[2]
    tt <- fastT(cbind(a.samples,b.samples),1:ai,(ai+1):(ai+bi),var.equal=FALSE)
    tt <- 2* pt(-abs(tt$z),df=ai+bi-2)
  }
  else {
    i <- 1:(dim(a.samples)[1])
    tt <- sapply(i,function(y) { t.test(a.samples[y,],b.samples[y,],paired=T)$p.val})
  }

  names(fc)       <- rownames(a.samples);
  rownames(means) <- rownames(a.samples);
  names(tt)       <- rownames(a.samples);


  pd <- pData(x);
  grp <- pd[,colnames(pd) == group];
  these.pd <- pd[is.element(grp,members),]

  
  return(new("PairComp",fc=fc,tt=tt,means=means,group=group,members=members,pData=these.pd,calculated.from=x))
}


"pairwise.comparison" <- function(x,group,members=NULL,spots=NULL,a.order=NULL,b.order=NULL,method="unlogged",logged=TRUE) {
  if(is.null(members)) {
    pd <- unique(as.character(pData(x)[,group]))
    if(is.null(pd)) {
      stop(paste("Can't find a group called",group));
    }      
    if(length(pd) != 2)  {
      stop("There must be exactly two groups for a pairwise comparison. Please specify which groups you want to compare.");
    }
    members <- pd;
  }
  if(!is.null(spots)) {
    pd <- pData(x);
    grp <- pd[,colnames(pd) == group];
    my.spots <- spots[,is.element(grp,members)]

    pmac    <- detection.p.val(my.spots);
    results <- get.fold.change.and.t.test(x,group,members,logged=logged,a.order=a.order,b.order=b.order,method=method);

    results@calls <- pmac$call;
  }
  else {
    results <- get.fold.change.and.t.test(x,group,members,logged=logged,a.order=a.order,b.order=b.order,method=method);
  }
  return(results); 
}




all.present.in.group <- function(x,group,members=NULL,calls,no="all") {
  
  pd      <- pData(x)
  grp     <- as.factor(pd[,colnames(pd) == group])
  if(is.null(members)) {
    members <- levels(grp)
  }
  colnames(calls) <- grp
  result  <- rep(FALSE,dim(calls)[1])
  for(i in members) {
    print(paste("Checking member",i,"in group: '",group,"'"))
    this.group <- calls[,grp %in% i]
    this.group <- array(sapply(this.group,
                        function(val) { if(val == "P") { 1 } else { 0 } } 
                       ),dim=dim(this.group));
    sums <- apply(this.group,1,sum)

    if(no=="all") { expected.number <- dim(this.group)[2] }
    else { expected.number <- no}
    ok   <- (sums >= expected.number) 
    result <- result | ok;
  }
  return(result)  
}



all.present <- function(x,calls,no = "all") {
  pd      <- pData(x)
    calls <- array(sapply(calls,
                   function(val) { if(val == "P") { 1 } else { 0 } } 
                  ),dim=dim(calls));
    sums <- apply(calls,1,sum)
    if(no=="all") { expected.number <- dim(calls)[2] }
    else { expected.number <- no}
    ok   <- (sums >= expected.number)
  return(ok)  
}


pairwise.filter <- function(object,min.exp=log2(100),min.exp.no=0, min.present.no=0,present.by.group=T,fc=1.0,tt=0.001) {
  
  if(class(object) != "PairComp") { stop("Can only filter an object of class 'PairComp'"); }
  x <- calculated.from(object)
  pass.fc              <- (abs(fc(object)) > fc);
  pass.tt              <- (tt(object) < tt );

  samples <- exprs(get.array.subset(x,group(object),members(object)));
  samples <- samples[names(fc(object)),]
  no.chips             <- length(colnames(samples));

  intensity.thresh     <- array(sapply(samples,function(x) { if(x > min.exp) { 1 } else { 0 } } ),dim=dim(samples));

  min.intensity.thresh <- rowSums(intensity.thresh)

  
  pass.intensity <- (min.intensity.thresh >= min.exp.no);

  if(nrow(calls(object))>0) {
    if(present.by.group) {
      pass.present <- all.present.in.group(object,group(object),members(object),calls(object),min.present.no);
    }
    else {
      pass.present <- all.present(object,calls(object),min.present.no);
    }
    return(object[(pass.fc & pass.tt & pass.intensity & pass.present),]);
  }
  else {
    return(object[(pass.fc & pass.tt & pass.intensity),]);
  }
}


plot.pairwise.comparison <- function(x,y=NULL,labels=colnames(means(x)),showPMA=TRUE,type="scatter",...) {
  if(type=="scatter") { .pcscatterplot(x,y,labels,showPMA,...) }
  else if(type=="volcano") { .volcanoplot(x,y,labels,showPMA,...) }
  else if(type=="ma") { .maplot(x,y,labels,showPMA,...) }
}


.pcscatterplot <- function(x,y,labels,showPMA,...) {
	pd <- pData(x)
	gp <- pd[,colnames(pd) == group(x)]
	me <- members(x)
	calls <- calls(x)
        d <- dim(calls)
	if(!(d[1] == 0 & d[2] == 0)) {
	  	callsA <- calls[,gp == me[1]]
	  	callsB <- calls[,gp == me[2]]
	  	sumsA <- apply(callsA,1,function(row) { 
                    apply(sapply(row,function(val) {
                         if(val == "P") {return(c(1,0,0))}
                         if(val == "M") {return(c(0,1,0))}
                         if(val == "A") {return(c(0,0,1))}
                    }),1,sum)
              })
	  	sumsB <- apply(callsB,1,function(row) { 
                    apply(sapply(row,function(val) {
                         if(val == "P") {return(c(1,0,0))}
                         if(val == "M") {return(c(0,1,0))}
                         if(val == "A") {return(c(0,0,1))}
                    }),1,sum)
              })
 	   	inA <- length(colnames(callsA))
 	   	inB <- length(colnames(callsB))
	  	allPA <- sumsA[1,] == inA
	  	allPB <- sumsB[1,] == inB
	  	allAA <- sumsA[3,] == inA
	  	allAB <- sumsB[3,] == inB
	  	allP  <- allPA & allPB
	  	allA  <- allAA & allAB
              AorB  <- (allPA | allPB) & !(allA | allP)
	  	
	  	trad.scatter.plot(means(x)[,1],means(x)[,2],xlab=labels[1],ylab=labels[2],col="yellow",...);
	  	trad.scatter.plot(means(x)[AorB,1],means(x)[AorB,2],add=T,col="orange");
	  	trad.scatter.plot(means(x)[allP,1],means(x)[allP,2],add=T,col="red");
	}
	else {
	  	trad.scatter.plot(means(x)[,1],means(x)[,2],xlab=labels[1],ylab=labels[2],col="light grey",...); 
 	}
	if(!missing(y)) {
    	  trad.scatter.plot(means(y)[,1],means(y)[,2],add=T,pch=1,col="blue");
	}

}


setMethod("plot",signature(x="PairComp",y="missing"),function(x,y,...) plot.pairwise.comparison(x,y,...))
setMethod("plot",signature(x="PairComp",y="PairComp"),function(x,y,...) plot.pairwise.comparison(x,y,...))


.volcanoplot <- function(x,y,labels,...) {
	pd <- pData(x)
	gp <- pd[,colnames(pd) == group(x)]
	me <- members(x)
	calls <- calls(x)
        d <- dim(calls)
	if(!(d[1] == 0 & d[2] == 0)) {
	  	callsA <- calls[,gp == me[1]]
	  	callsB <- calls[,gp == me[2]]
	  	sumsA <- apply(callsA,1,function(row) { 
                    apply(sapply(row,function(val) {
                         if(val == "P") {return(c(1,0,0))}
                         if(val == "M") {return(c(0,1,0))}
                         if(val == "A") {return(c(0,0,1))}
                    }),1,sum)
              })
	  	sumsB <- apply(callsB,1,function(row) { 
                    apply(sapply(row,function(val) {
                         if(val == "P") {return(c(1,0,0))}
                         if(val == "M") {return(c(0,1,0))}
                         if(val == "A") {return(c(0,0,1))}
                    }),1,sum)
              })
 	   	inA <- length(colnames(callsA))
 	   	inB <- length(colnames(callsB))
	  	allPA <- sumsA[1,] == inA
	  	allPB <- sumsB[1,] == inB
	  	allAA <- sumsA[3,] == inA
	  	allAB <- sumsB[3,] == inB
	  	allP  <- allPA & allPB
	  	allA  <- allAA & allAB
                AorB  <- (allPA | allPB) & !(allA | allP)

	   plot(fc(x),log2(tt(x)),xlab=paste(labels[2],"<-- fold change -->",labels[1]),ylab="p-score",pch=20,col="yellow")
	   points(fc(x)[AorB],log2(tt(x)[AorB]),pch=20,col="orange")
   	   points(fc(x)[allP],log2(tt(x)[allP]),pch=20,col="red")
       }
       else {
	   plot(fc(x),log2(tt(x)),xlab=paste(labels[2],"<-- fold change -->",labels[1]),ylab="p-score",pch=20,col="light grey")
       }
   if(!missing(y)) {
     points(fc(y),log2(tt(y)),col="blue",pch=1);
   }
}


.maplot<- function(x,y,labels,...) {
	m <- fc(x);
        a <- (means(x)[,1] + means(x)[,2])/2
	pd <- pData(x)
	gp <- pd[,colnames(pd) == group(x)]
	me <- members(x)
	calls <- calls(x)
        d <- dim(calls)
	if(!(d[1] == 0 & d[2] == 0)) {
	  	callsA <- calls[,gp == me[1]]
	  	callsB <- calls[,gp == me[2]]
	  	sumsA <- apply(callsA,1,function(row) { 
                    apply(sapply(row,function(val) {
                         if(val == "P") {return(c(1,0,0))}
                         if(val == "M") {return(c(0,1,0))}
                         if(val == "A") {return(c(0,0,1))}
                    }),1,sum)
              })
	  	sumsB <- apply(callsB,1,function(row) { 
                    apply(sapply(row,function(val) {
                         if(val == "P") {return(c(1,0,0))}
                         if(val == "M") {return(c(0,1,0))}
                         if(val == "A") {return(c(0,0,1))}
                    }),1,sum)
              })
 	   	inA <- length(colnames(callsA))
 	   	inB <- length(colnames(callsB))
	  	allPA <- sumsA[1,] == inA
	  	allPB <- sumsB[1,] == inB
	  	allAA <- sumsA[3,] == inA
	  	allAB <- sumsB[3,] == inB
	  	allP  <- allPA & allPB
	  	allA  <- allAA & allAB
                AorB  <- (allPA | allPB) & !(allA | allP)

	   plot(a,m,xlab="A",ylab="M",pch=20,col="yellow")
	   points(a[AorB],m[AorB],pch=20,col="orange")
   	   points(a[allP],m[allP],pch=20,col="red")
       }
       else {
	   plot(a,m,xlab="A",ylab="M",pch=20,col="light grey")
       }
   if(!missing(y)) {
	m <- fc(y);
        a <- (means(y)[,1] + means(y)[,2])/2

     points(a,m,col="blue",pch=1);
   }

}


