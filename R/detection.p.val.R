detection.p.val <- function(x, tau = NULL,calls=TRUE,alpha1=NULL,alpha2=NULL,ignore.saturated=TRUE) {
  if(is.null(tau)) {
    if(!qc.ok()) {
      setQCEnvironment(cleancdfname(cdfName(x)))
    }
    tau <- qc.get.tau()
  }
  if(is.null(alpha1)) {
    if(!qc.ok()) {
      setQCEnvironment(cleancdfname(cdfName(x)))
    }
    alpha1 <- qc.get.alpha1()
  }
  if(is.null(alpha2)) {
    if(!qc.ok()) {
      setQCEnvironment(cleancdfname(cdfName(x)))
    }
    alpha2 <- qc.get.alpha2()
  }
  if(class(x) != "AffyBatch") { stop("detection.p.val() should be called on an AffyBatch object.\nSee ?detection.p.val for more details.") }
  if(alpha1 < 0)      {stop("alpha1 must be  > 0 "); }
  if(alpha1 > alpha2) {stop("alpha2 must be  > alpha1 "); }
  if(alpha2 > 1)      {stop("alpha2 must be  <1 "); }

  pms <-as.matrix(pm(x));
  mms <-as.matrix(mm(x));

  # Saturation:
  # shouldn't be a problem with new scanners or those that have had an engineer visit
  if(ignore.saturated) { sat <- 46000; }
  else { sat <- -1; }
  
  pns <- probeNames(x);
  o <- order(pns)
  pns <- pns[o]
  pms <- pms[o,,drop=FALSE]
  mms <- mms[o,,drop=FALSE]
  unique.pns <- sort(unique(pns));

  p<-sapply(1:length(pms[1,]),function(x) { 
    .C("DetectionPValue",as.double(pms[,x]),as.double(mms[,x]),as.character(pns),as.integer(length(mms[,x])),
	as.double(tau),as.double(sat),dpval=double(length(unique.pns)),length(unique.pns),PACKAGE="simpleaffy")$dpval;
  });
  rownames(p) <- unique.pns;
  colnames(p) <- paste(sampleNames(x),".detection.p.val",sep="");
  if(!calls) { 
    l <- list(detection.p.values=p); 
  }
  else       {
    calls <- sapply(p,function(y) { if(y < alpha1) { return("P") } else { if(y < alpha2) { return("M") } else { return("A") }}});
    calls <- matrix(calls,nrow=nrow(p),ncol=ncol(p));
    colnames(calls) <- paste(sampleNames(x),".present",sep="");
    rownames(calls) <- rownames(p)
    l<- list(pval=p,call=calls);
    return(l); 
  }

}     

