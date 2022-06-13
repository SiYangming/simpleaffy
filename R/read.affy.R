read.affy <- function(covdesc="covdesc",path=".",...) {
  samples <- read.AnnotatedDataFrame( paste(path,covdesc,sep="/"),sep="");
  files.to.read <- rownames(pData(samples));
  files.to.read <- paste(path,files.to.read,sep="/")
  eset <- ReadAffy(filenames=files.to.read,...);
  newPhenoData <- cbind(pData(eset),pData(samples)[rownames(pData(eset)),]);
  colnames(newPhenoData) <- c(colnames(pData(eset)),colnames(pData(samples)));
  tmp <- as.list(colnames(newPhenoData));
  names(tmp) <- colnames(newPhenoData);
  newPhenoData <- as(newPhenoData,"AnnotatedDataFrame")

  phenoData(eset) <- newPhenoData;	
  return(eset);
}


"read.affy.mixed" <-
function(covdesc="covdesc",path=".",...) {
  samples <- read.AnnotatedDataFrame( paste(path,covdesc,sep="/"),sep="");

  files.to.read <- rownames(pData(samples));
  files.to.read <- paste(path,files.to.read,sep="/")
  esets <- lapply(files.to.read,function(x) {ReadAffy(filenames=x,...)})

  nms   <- sapply(esets,cdfName)
  u.nms <- unique(nms)
  merged <- list()
  for(i in 1:length(nms)) {
   so.far <- merged[nms[i]];
   if(is.null(so.far[[1]])) {
     merged[[nms[i]]]<-esets[[i]]

   }
   else {
     merged[[nms[i]]] <- merge.AffyBatch(so.far[[1]],esets[[i]])
   }
  }

  nms <- unique(sapply(merged,cdfName))

  acc <- list(probeNames(merged[[1]]))
  nm  <- cdfName(merged[[1]])[1]

  idxs <- list()
  idxs[nm] <- acc


  common <- as.vector(unlist(acc))
  for(i in 2:length(merged)) {
  
    nm  <- cdfName(merged[[i]])[1]

    if(is.null(idxs[nm][[1]])) {
      acc <- list(probeNames(merged[[i]]))

      idxs[nm] <- acc
      common <- intersect(common,as.vector(unlist(acc)))
    }
  }

  smallest <- 1;

  for(i in 2:length(merged)) {
    if(length(probeNames(merged[[i]])) < length(probeNames(merged[[smallest]]))) { smalleset <- i }
  }
  
  all.smallest <- unique(probeNames(merged[[smallest]]))

  template <- merged[[smallest]][1]
  not.shared <- setdiff(all.smallest,common)
  idxs <- indexProbes(template,"both")
  e <- exprs(template)
  e[unlist(idxs[is.element(all.smallest,not.shared)]),] <- 0


  exprs(template) <- e

  result <- template;
  for(i in 2:length(esets)) {
    result <- merge.AffyBatch(result,template)
  }

  to.change <- exprs(result)
 
  idxs.small <- unlist(idxs[common])
  o <- order(names(idxs.small))
  idxs.small <- idxs.small[o];
  
  changing <- 1;

  for(i in 1:length(merged)) {
    idxs.large <- unlist(indexProbes(merged[[i]],"both")[common])
    o <- order(names(idxs.large))
    idxs.large <- idxs.large[o];

    for(j in 1:length(merged[[i]])) {
	  to.copy <- exprs(merged[[i]][j])
	  to.change[idxs.small,changing] <- to.copy[idxs.large,1]
          changing <- changing + 1;
    } 
  }
  exprs(result) <- to.change;
  rn <- paste(path,rownames(pData(samples)),sep="/")
  rownames(pData(samples)) <- rn 
  newPhenoData <- pData(samples)[rn,];
  tmp <- as.list(colnames(newPhenoData));
  names(tmp) <- colnames(newPhenoData);
  newPhenoData <- new("AnnotatedDataFrame")
  pData(newPhenoData) <- newPhenoData
  varLabels(newPhenoData) <- tmp

  result@phenoData <- newPhenoData;
  result@notes<-"This expression set was produced by merging different chip types. Be careful!"	
  cat("done!\n Be careful how you use the results from this function...\n");
  return(result)
}
