# Provides a wrapper around some of the various expression calling
# algorithms on offer. 



"call.exprs" <-
function(x,algorithm="rma",do.log = TRUE,sc=100,method=NA) {

  if(class(x) != "AffyBatch") { 
	stop("call.exprs() should be called on an an AffyBatch object.\nSee ?call.exprs for more detail.");
  }

  if(algorithm == "rma-R") {       # use RMA
    if(is.na(method)) {
      method<-"quantiles";
    }
    tmp <- expresso(x, bgcorrect.method="rma",
                      normalize.method=method,pmcorrect.method="pmonly",
                      summary.method="avgdiff");
    if(!do.log) {
      exprs(tmp) <- log2(exprs(tmp));
    }
    return(tmp);
  }
   if(algorithm == "rma") {       # use RMA
    tmp <- rma(x);
    if(!do.log) {
      exprs(tmp) <- 2^exprs(tmp);
    }
    return(tmp);
  }
  else if(algorithm == "gcrma") {       # use GCRMA
    if(is.na(method)) {
      method<-"quantiles";
    }
    tmp <- gcrma(x);
    if(!do.log) {
      exprs(tmp) <- 2^exprs(tmp);
    }
    return(tmp);
  }

  else if(algorithm == "mas5-R") { # use expresso MAS5.0
     if(is.na(method)) {
     tmp1 <- expresso(x, normalize=FALSE,bgcorrect.method="mas",
                        pmcorrect.method="mas",summary.method="mas");
        tmp  <- affy.scalevalue.exprSet(tmp1,sc=sc);
     }
     else {
     tmp1 <- expresso(x, normalize.method=method,bgcorrect.method="mas",
                        pmcorrect.method="mas",summary.method="mas");
        tmp  <- affy.scalevalue.exprSet(tmp1,sc=sc);
     }
     if(do.log) {
       exprs(tmp) <- log2(exprs(tmp));
       preproc(tmp)$sfs=apply(2^(exprs(tmp) - log2(exprs(tmp1))),2,mean)
       preproc(tmp)$tgt=sc
     }
     else {
       preproc(tmp)$sfs=apply(2^(exprs(tmp) - log2(exprs(tmp1))),2,mean)
       preproc(tmp)$tgt=sc

     }
     return(tmp);
  }
  else if(algorithm == "mas5") { # use Simpleaffy MAS5.0
    tmp <- justMAS(x,tgt=sc);
    if(!do.log) {
      exprs(tmp) <- 2^exprs(tmp);
    }
    return(tmp);
  }
  else {
    stop(paste("Don't know how to compute algorithm '",algorithm,"'.\nSee ?call.exprs() for more details."));
  }
}
