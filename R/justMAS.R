bg.correct.sa <- function(unnormalised, grid=c(4,4)) {
res         <- unnormalised;
pms         <- unique(unlist(pmindex(res))) - 1 # C counts from 0
mms         <- unique(unlist(mmindex(res))) - 1 # C counts from 0
all         <- c(pms,mms)
intensities <- intensity(res)
rws <- nrow(res)
cls <- ncol(res)
for(no in 1:length(res)){
  this.array <- intensities[,no];
  result <- .C("bgmas",as.integer(as.vector(all)),as.integer(length(all)),
       as.double(as.vector(this.array)),as.integer(length(this.array)),
       as.integer(rws),
       as.integer(cls),
       as.integer(grid[1]),as.integer(grid[2]),
       zonebg=double(grid[1] * grid[2]),
       zonesd=double(grid[1] * grid[2]),corrected=double(length(this.array)),PACKAGE="simpleaffy");
       intensities[,no] <- result$corrected;
  }
  intensity(res) <- intensities;
  return(res);
}

justMAS     <- function(unnormalised,tgt=100,scale=TRUE) {
  ct <- 0.03;
  st <- 10.0;
########################### BACKGROUND
  bgc <- bg.correct.sa(unnormalised);

  pms <-as.matrix(pm(bgc))
  mms <-as.matrix(mm(bgc))
  pns <- probeNames(bgc);
  o <- order(pns)
  pns <- pns[o]
  pms <- pms[o,,drop=FALSE]
  mms <- mms[o,,drop=FALSE]
  unique.pns <- sort(unique(pns));

########################### SUMMARIES
  expression.calls<-sapply(1:length(pms[1,]),function(x) { 
    .C("GetExpressionLevels",as.double(pms[,x]),as.double(mms[,x]),as.character(pns),as.integer(length(mms[,x])),
	as.double(ct),as.double(st),exprs=double(length(unique.pns)),length(unique.pns),PACKAGE="simpleaffy")$exprs;
  });
  rownames(expression.calls) <- unique.pns;
  colnames(expression.calls) <- paste(sampleNames(bgc))

  res <- new("ExpressionSet", 
              exprs       = expression.calls,
              phenoData   = as(phenoData(bgc),"AnnotatedDataFrame"),
              annotation  = annotation(bgc)) 
  ##FIXME: bug in ExpressionSet initialization code
              #notes       = notes(bgc));

########################### SCALING
  if(scale) {
   sfs <- double(length(expression.calls[1,]));


   for(x in 1:length(expression.calls[1,])) {
     vals <- sort(2^expression.calls[,x]);
     l <- length(vals);
     frm <- 0.02 *l;
     to  <- 0.98 *l;
     sf  <- tgt/mean(vals[frm:to]);
     expression.calls[,x] <- log2((2^expression.calls[,x]) * sf)
     sfs[x] <- sf; 
   }
   exprs(res) <- expression.calls;
   expData <- experimentData(res)
   expData@preprocessing <- c(expData@preprocessing,list(sfs = unlist(sfs),tgt = tgt))
   experimentData(res) <- expData
  }
  else {
   expData <- experimentData(res)
   expData@preprocessing <- c(expData@preprocessing,list(sfs = stop("Arrays were not scaled"), stop("Arrays were not scaled")))
   experimentData(res) <- expData
  }

  return(res);
}


.mas5     <- function(unnormalised,normalize=TRUE,sc=500,analysis="absolute") {
  ct <- 0.03;
  st <- 10.0;
########################### BACKGROUND
  if(normalize) {
    bgc <- bg.correct.sa(unnormalised);
  }
  pms <-as.matrix(pm(bgc))
  mms <-as.matrix(mm(bgc))
  pns <- probeNames(bgc);
  o <- order(pns)
  pns <- pns[o]
  pms <- pms[o,,drop=FALSE]
  mms <- mms[o,,drop=FALSE]
  unique.pns <- sort(unique(pns));
########################### SUMMARIES
  expression.calls<-sapply(1:length(pms[1,]),function(x) { 
    .C("GetExpressionLevels",as.double(pms[,x]),as.double(mms[,x]),as.character(pns),as.integer(length(mms[,x])),
	as.double(ct),as.double(st),exprs=double(length(unique.pns)),length(unique.pns),PACKAGE="simpleaffy")$exprs;
  });
  rownames(expression.calls) <- unique.pns;
  colnames(expression.calls) <- paste(sampleNames(bgc))
########################### SCALING
  sfs <- double(length(expression.calls[1,]));
  for(x in 1:length(expression.calls[1,])) {
    vals <- sort(2^expression.calls[,x]);
    l <- length(vals);
    frm <- 0.02 *l;
    to  <- 0.98 *l;
    sf  <- sc/mean(vals[frm:to]);
    expression.calls[,x] <- log2((2^expression.calls[,x]) * sf)
    sfs[x] <- sf; 
  }
  res <- new("ExpressionSet", 
             exprs       = expression.calls,
             phenoData   = as(phenoData(bgc),"AnnotatedDataFrame"),
             annotation  = annotation(bgc))
   ##FIXME: come back and put this back after ExpressionSet bug is fixed
   ##          notes       = notes(bgc));
   expData <- experimentData(res)
   expData@preprocessing <- c(expData@preprocessing,list(sfs = stop("Arrays were not scaled"), stop("Arrays were not scaled")))
   experimentData(res) <- expData
}

