#holds the results of a QC comparison
setClass("QCStats",representation(scale.factors="numeric",target="numeric",percent.present="numeric",average.background="numeric",minimum.background="numeric",maximum.background="numeric",spikes="matrix",qc.probes="matrix",bioBCalls="character",arraytype="character"));

#accessor methods
setGeneric("sfs", function(object) standardGeneric("sfs"))
setMethod("sfs","QCStats",function(object) object@scale.factors)

setGeneric("target", function(object) standardGeneric("target"))
setMethod("target","QCStats",function(object) object@target)

setGeneric("percent.present", function(object) standardGeneric("percent.present"))
setMethod("percent.present","QCStats",function(object) object@percent.present)

setGeneric("avbg", function(object) standardGeneric("avbg"))
setMethod("avbg","QCStats",function(object) object@average.background)

setGeneric("minbg", function(object) standardGeneric("minbg"))
setMethod("minbg","QCStats",function(object) object@minimum.background)

setGeneric("maxbg", function(object) standardGeneric("maxbg"))
setMethod("maxbg","QCStats",function(object) object@maximum.background)

setGeneric("spikeInProbes", function(object) standardGeneric("spikeInProbes"))
setMethod("spikeInProbes","QCStats",function(object) object@spikes)

setGeneric("qcProbes", function(object) standardGeneric("qcProbes"))
setMethod("qcProbes","QCStats",function(object) object@qc.probes)

setGeneric("arrayType", function(object) standardGeneric("arrayType"))
setMethod("arrayType","QCStats",function(object) object@arraytype)


# alter the QC environment - effectively accessor functions into the environment

.qc.is.empty <- function() {
  get("empty",.qcEnv)
}

.qc.set.empty <- function(empty) {
  assign("empty",empty,envir=.qcEnv)
}


qc.set.array <- function(name) {
  .qc.set.empty(FALSE)
  assign("array",name,envir=.qcEnv)
}

qc.get.array <- function() {
  .qc.test()
  get("array",.qcEnv)
}

qc.get.tau <- function() {
  .qc.test()
  0.015;
}

qc.set.alpha1 <- function(value) {
  .qc.set.empty(FALSE)
  assign("alpha1",value,envir=.qcEnv)
}

qc.get.alpha1 <- function() {
  .qc.test()
  get("alpha1",.qcEnv)
}

qc.set.alpha2 <- function(value) {
  .qc.set.empty(FALSE)
  assign("alpha2",value,envir=.qcEnv)
}

qc.get.alpha2 <- function() {
  .qc.test()
  get("alpha2",.qcEnv)
}

qc.get.spikes <- function() {
  .qc.test()
  get("spikes",.qcEnv)
}

qc.get.spike <- function(name) {
  .qc.test()
  spikes <- get("spikes",.qcEnv)
  spikes[[name]]
}

qc.add.spike <- function(name,probeset) {
  .qc.set.empty(FALSE)
  spikes <- qc.get.spikes()
  names  <- names(spikes)
  if(name %in% names) {
    spikes[name] <- probeset
  }
  else {
    spikes <- c(spikes,probeset)
    names  <- c(names,name)
    names(spikes) <- names
  }
  assign("spikes",spikes,envir=.qcEnv)
}

qc.get.probes <- function() {
  .qc.test()
  get("probes",.qcEnv)
}

qc.get.probe <- function(name) {
  .qc.test()
  probes <- get("probes",.qcEnv)
  probes[[name]]
}

qc.add.probe <- function(name,probeset) {
  .qc.set.empty(FALSE)
  probes <- qc.get.probes()
  names  <- names(probes)
  if(name %in% names) {
   probes[name] <- probeset
  }
  else {
    probes <- c(probes,probeset)
    names  <- c(names,name)
    names(probes) <- names
  }
  assign("probes",probes,envir=.qcEnv)
}

qc.get.ratios <- function() {
  get("ratios",.qcEnv)
}

qc.get.ratio <- function(name) {
  .qc.test()
  probes <- get("ratios",.qcEnv)
  probes[[name]]
}

qc.add.ratio <- function(name,probeset1,probeset2) {
  .qc.set.empty(FALSE)
  ratios <- qc.get.ratios()
  names  <- names(ratios)
  if(name %in% names) {
    ratios[[name]] <- c(probeset1,probeset2)
  }
  else {
    if(length(ratios) == 0) {
      ratios <- list(c(probeset1,probeset2))
    }
    else {
      ratios <- c(ratios,list(c(probeset1,probeset2)))
    }
    names  <- c(names,name)
    names(ratios) <- names
  }
  assign("ratios",ratios,envir=.qcEnv)
}

qc.ok <- function() {
  !(.qc.is.empty())
}

.qc.test <- function() {
  if(.qc.is.empty()) {
    stop("QC environment is empty.\n See 'setQCEnvironment' for more details.\n")
  }
}
    
# read a file and use this to set up the QC Environment
# Expects a (whitespace delimited) file set out like this:
#
#
#    array <arrayname>
#    ratio <rationame> <pset1> <pset2>
#    ratio <rationame> <pset1> <pset2>
#    ratio <rationame> <pset1> <pset2>
#    ...
#    spk <spikename> <pset>
#    spk <spikename> <pset>
#    spk <spikename> <pset>
#    spk <spikename> <pset>
#    ...
#    alpha1 <value>
#    alpha2 <value>
#
# where spikename can be one of bioB, bioC, bioD or creX

qc.read.file <- function(fn) {
  .createEmptyQCEnvironment()
  fl <- file(fn,"r")
  if(!isOpen(fl)) { stop(paste("Couldn't open file",fn,"\n.")) }
  lines <- readLines(fl)
  lines <- strsplit(lines,"\\s+")
  for(l in lines) {
    a <- NULL
    if(length(l) > 0) {
      switch(l[1],
             array  = .qc.file.setQCArray(l),
             spk    = .qc.file.addQCSpike(l,a),
             ratio  = .qc.file.addQCRatio(l,a),
             alpha1 = .qc.file.setAlpha1(l,a),
             alpha2 = .qc.file.setAlpha2(l,a))
    }
  }
  close(fl)
}

.qc.file.setQCArray <- function(toks) {
  if(length(toks) != 2) {
    stop("Array name should be a single string.");
  }
  qc.set.array(toks[2])
}

#Adds or replaces the probeset for the specified QC spike on the array, 'a'
.qc.file.addQCSpike <- function(toks,a) {
  if(length(toks) != 3) {
    stop(paste("Error parsing file\nExpecting: spk '[bioB|bioC|bioD|creX] <probesetid>'\nGot: '",paste(toks,collapse=" "),"'.\n"));
  }
  if((toks[2] != "bioB") &
     (toks[2] != "bioC") &
     (toks[2] != "bioD") &
     (toks[2] != "creX"))  {
    stop(paste("Error parsing file\nSpike name must be one of  'bioB', 'bioC', 'bioD' or 'creX'\nGot: '",paste(toks,collapse=" "),"'.\n"));
  }
  qc.add.spike(toks[2],toks[3])
}

#store in the ratios list in .qcEnv each pairwise comparison. Add the relevant probesets to qc.probes, if they are not already there
.qc.file.addQCRatio <- function(toks,a) {
  if(length(toks) != 4) {
    stop(paste("Error parsing file\nExpecting: ratio <probeset1name>/<probeset2name> <probesetid1> <probesetid2>'\nGot: ",c(toks),"'.\n"));
  }
  name <- toks[2]
  ps1 <- toks[3]
  ps2 <- toks[4]
  probenames <- strsplit(name,"/")[[1]]
  if(length(probenames)!=2) {
    stop(paste("Error parsing file\nExpecting rationame: '<probeset1name>/<probeset2name>'\nGot: ",name,"'.\n"));
  }
  qc.add.probe(probenames[1],ps1)
  qc.add.probe(probenames[2],ps2)
  qc.add.ratio(name,ps1,ps2)
}

.qc.file.setAlpha1 <- function(toks,a) {
  if(length(toks) != 2) {
    stop(paste("Error parsing file\nExpecting: alpha1 <value>'\nGot: ",paste(toks),"'.\n"));
  }
  v <- as.numeric(toks[2])
  if(is.na(v)) {
    stop(paste("Error parsing file\nalpha1 value must be a number'\nGot: ",paste(toks),"'.\n"));
  }

  qc.set.alpha1(v)
}


.qc.file.setAlpha2 <- function(toks,a) {
  if(length(toks) != 2) {
    stop(paste("Error parsing file\nExpecting: alpha2 <value>'\nGot: ",paste(toks),"'.\n"));
  }
  v <- as.numeric(toks[2])
  if(is.na(v)) {
    stop(paste("Error parsing file\nalpha1 value must be a number'\nGot: ",paste(toks),"'.\n"));
  }

  qc.set.alpha2(v)
}

.qcEnv <- new.env(parent=emptyenv())

#initialize the QC environment 
.createEmptyQCEnvironment <- function() {
  assign("empty",TRUE,envir=.qcEnv)
  assign("array",NULL,envir=.qcEnv)
  assign("alpha1",NULL,envir=.qcEnv)
  assign("alpha2",NULL,envir=.qcEnv)
  assign("spikes",vector(),envir=.qcEnv)
  assign("probes",vector(),envir=.qcEnv)
  assign("ratios",list(),envir=.qcEnv)
}


qc.have.params <- function(name) {
  fn <- .get.array.def.file(name)
  !(fn == "")
}


.initializeQCEnvironment <- function() {
    .createEmptyQCEnvironment();
}

.initializeQCEnvironment();


# set up the QC environment for a given array type
# if path specified, look there. Otherwise
# look in extdata to see if there's a file there
setQCEnvironment <- function(array,path=NULL) {
  if(is.null(path)) {
    fn <- .get.array.def.file(array)
    if(fn != "") {
      qc.read.file(fn)
    }
    else {
      stop(paste("Could not find array definition file '",paste(array,"qcdef",sep="."),"'. Simpleaffy does not know the QC parameters for this array type.\nSee the package vignette for details about how to specify QC parameters manually.\n"))
    }
  }
  else {
    fn <- file.path(path,paste(array,"qcdef",sep="."))
    if(fn != "") {
      qc.read.file(fn)
    }
    else {
      stop(paste("Could not find array definition file: '",paste(array,"qcdef",sep="."),"' in the directory specified:",path,"\n"))
    }
  }
}

.get.array.def.file <- function(array) {
  fn <- paste(array,"qcdef",sep=".")
  system.file("extdata",fn,package="simpleaffy")
}



.getRatios <- function(x) {
   vals <- qcProbes(x);
   if(!qc.ok()) {
     setQCEnvironment(arrayType(x))
   }
   to.calculate <- qc.get.ratios();
   r <- sapply(to.calculate,function(a) { vals[,a[1]] - vals[,a[2]]})
   return(r)
}

setGeneric("ratios", function(object) standardGeneric("ratios"))
setMethod("ratios","QCStats",function(object) .getRatios(object))


qc.affy <- function(unnormalised,normalised=NULL,tau=0.015,logged=TRUE,cdfn=cdfName(unnormalised)) {
   verbose <- getOption("verbose")
   if(.qc.is.empty()) {
     cdfn <- cleancdfname(cdfn)


     if(verbose){cat(paste("Looking cdf file for:",cdfn,"\n"))}

     setQCEnvironment(cdfn)
  }
  else {
   if(qc.get.array() !=  cleancdfname(cdfn)) {
	warning(paste("CDF Environment name '", qc.get.array(), "' does not match cdfname '", cleancdfname(cdfn),"'"))
   }
  }
  if(is.null(normalised)) {
    if(verbose){cat(paste("Preprocessing expression data using mas5\n"))}
    normalised <- call.exprs(unnormalised,"mas5");
  }

  x <- exprs(normalised);

  det <- detection.p.val(unnormalised,tau=tau,alpha1=qc.get.alpha1(),alpha2=qc.get.alpha2());

   #change
  dpv<-apply(det$call,2,function(x) {
            x[x!="P"] <- 0;
	    x[x=="P"] <- 1;
            x<-as.numeric(x);			       
            return(100 * sum(x)/length(x));
  });

  sfs    <- experimentData(normalised)@preprocessing$sfs;
  target <- experimentData(normalised)@preprocessing$tgt;

  if(!logged) { x <- log2(x); }

  bgsts<-.bg.stats(unnormalised)$zonebg

  meanbg <- apply(bgsts,1,mean);

  minbg  <- apply(bgsts,1,min);

  maxbg  <- apply(bgsts,1,max);

  stdvbg <- sqrt(apply(bgsts,1,var));




  #get the probenames for the QC probes for this chip
  qc.probenames <- qc.get.probes();
  
  qc.probe.vals <- rbind(c(),(sapply(qc.probenames, function(y) {x[y,]})))
  rownames(qc.probe.vals) <- colnames(x);
  colnames(qc.probe.vals) <- qc.probenames
  spike.probenames <- qc.get.spikes();
  spike.vals <- rbind(c(),(sapply(spike.probenames, function(y) {x[y,]})))

  rownames(spike.vals) <- colnames(x);
  colnames(spike.vals) <- spike.probenames
   
  bb <- spike.probenames["bioB"]

  if(!is.na(bb)) {
    biobcalls <- det$call[bb,]
  }
  else {
     biobcalls <- NULL
  }

  return(new("QCStats",scale.factors=sfs,target=target,percent.present=dpv,average.background=meanbg,minimum.background=minbg,maximum.background=maxbg,
              spikes=spike.vals,qc.probes=qc.probe.vals,bioBCalls=biobcalls,arraytype=cdfn));
}


setGeneric("qc", function(unnormalised,...) standardGeneric("qc"))
setMethod("qc","AffyBatch",function(unnormalised,...) qc.affy(unnormalised,...)) 



.bg.stats <- function(unnormalised, grid=c(4,4)) {
pms         <- unlist(pmindex(unnormalised))
mms         <- unlist(mmindex(unnormalised))
all         <- c(pms,mms)
intensities <- exprs(unnormalised)
rws <- nrow(unnormalised)
cls <- ncol(unnormalised)
zonebg <- c();
zonesd <- c();
for(no in 1:length(unnormalised)){
  this.array <- intensities[,no];
  result <- .C("bgmas",as.integer(as.vector(all)),as.integer(length(all)),
       as.double(as.vector(this.array)),as.integer(length(this.array)),
       as.integer(rws),
       as.integer(cls),
       as.integer(grid[1]),as.integer(grid[2]),
       zonebg=double(grid[1] * grid[2]),
       zonesd=double(grid[1] * grid[2]),corrected=double(length(this.array)),PACKAGE="simpleaffy");
  zonesd <- rbind(zonesd, result$zonesd);
  zonebg <- rbind(zonebg, result$zonebg);
  }
  colnames(zonesd) <- paste("zone",1:16,sep=".");
  colnames(zonebg) <- paste("zone",1:16,sep=".");
  rownames(zonesd) <- sampleNames(unnormalised);
  rownames(zonebg) <- sampleNames(unnormalised);
  return(list(zonebg=zonebg,zonesd=zonesd))
}


plot.qc.stats<-function(x,fc.line.col="black",sf.ok.region="light blue",chip.label.col="black",sf.thresh = 3.0,gdh.thresh = 1.25,ba.thresh = 3.0,present.thresh=10,bg.thresh=20,label=NULL,main="QC Stats",usemid=F,spread=c(-8,8),cex=1,...) {
  old.par <- par()
  par(mai=c(0,0,0,0))
  sfs    <- log2(sfs(x))

  n      <- length(sfs)

  meansf <- mean(sfs)

  dpv <- percent.present(x)
  dpv <- (round(100*dpv))/100;

  abg <- avbg(x)
  abg <- (round(100*abg))/100;
	
  if(is.null(label)) { label <- names(maxbg(x)) }
  d1 <- 0.0;
  d2 <- 0.0;
  d3 <- 0.0;

  for(i in 1:n) {
    for(j in 1:n) { 
      d1 <- max(abs(sfs[i] - sfs[j]),d1);
      d2 <- max(abs(dpv[i] - dpv[j]),d2);
      d3 <- max(abs(abg[i] - abg[j]),d3);
    }
  }

  # set up plotting area - a column for array names next to a column for the QC

  m <- matrix(c(4,2,1,3) ,nrow=2,ncol=2)
  layout(m,c(1,2),c(0.1,1))
  # the title
  if(is.null(main)) { main="" }
  plot(0,0,xlim=range(0,1),ylim=range(0,1),type="n",yaxs="i",xaxt="n",yaxt="n",bty="n")
  text(0.5,0.5,labels=main,adj=0,cex=cex*2)

  # write out the array names

  plot(0,0,xlim=range(0,1),ylim=range(-1,n),type="n",yaxs="i",xaxt="n",yaxt="n",bty="n")
  text(1,(1:n)-0.5,labels=label,adj=1,cex=cex)
  plot(0,0,xlim=spread,ylim=c(-1,n),type="n",xaxs="i",yaxs="i",xaxt="n",yaxt="n",bty="n")

  x1 <- (sf.thresh/2.0 +  meansf)
  y1 <- 0
  x2 <- (-sf.thresh/2.0 +  meansf)
  y2 <- n

  polygon(c(x1,x2,x2,x1),c(y1,y1,y2,y2),col=sf.ok.region,border=sf.ok.region);
  lines(c(0,0),c(0,n),lty=1,col=fc.line.col)
  lines(c(-1,-1),c(0,n),lty=2,col="grey")
  lines(c(-2,-2),c(0,n),lty=2,col="grey")
  lines(c(-3,-3),c(0,n),lty=2,col=fc.line.col)
  lines(c(1,1),c(0,n),lty=2,col="grey")
  lines(c(2,2),c(0,n),lty=2,col="grey")
  lines(c(3,3),c(0,n),lty=2,col=fc.line.col)
  text(3,-1,"3",pos=3,col=fc.line.col,cex=cex)
  text(2,-1,"2",pos=3,col=fc.line.col,cex=cex)
  text(1,-1,"1",pos=3,col=fc.line.col,cex=cex)
  text(-3,-1,"-3",pos=3,col=fc.line.col,cex=cex)
  text(-2,-1,"-2",pos=3,col=fc.line.col,cex=cex)
  text(-1,-1,"-1",pos=3,col=fc.line.col,cex=cex)
  text(0,-1,"0",pos=3,col=fc.line.col,cex=cex)

  rats <- ratios(x);
  if(!usemid) {
    gdh <- rats[,3];
    ba  <- rats[,1];
  }
  else {
    gdh <- rats[,4];
    ba  <- rats[,2];
  }

  bb  <- x@bioBCalls

  for(i in 1:n) {
    x1<-spread[1]
    x2<-spread[2]
    y1<-i-1;
    y2<-i-1;
    lines(c(x1,x2),c(y1,y2),lty=2,col="light grey")
    if(d1 > sf.thresh) { col = "red" } else {col="blue"}
     x1 <- sfs[i]
     y1 <- i-0.25
     lines(c(0,x1),c(y1,y1),col=col);

     points(x1,y1,col=col,pch=20);
     x2 <- gdh[i]
     y2 <- i-0.5;
     if(gdh[i] > gdh.thresh) { col = "red" } else {col="blue"}	
     points(x2,y2,pch=1,col=col);

     x2 <- ba[i];
     y2 <- i-0.5;
     if(ba[i] > ba.thresh) { col = "red" } else {col="blue"}	
     points(x2,y2,pch=2,col=col);

     if(d2 > present.thresh) { col = "red" } else {col="blue"}
     x2 <- spread[1]
     y2 <- i-0.25
     dpvs<-paste(dpv[i],"%",sep="")
     text(x2,y2,label=dpvs,col=col,pos=4,cex=cex);
     if(d3 > bg.thresh) { col = "red" } else {col="blue"}
     x2 <- spread[1]
     y2 <- i-0.75
     text(x2,y2,label=abg[i],col=col,pos=4,cex=cex);
     if(bb[i]!="P") {
       text(0,i-1,label="bioB",col="red",cex=cex);
     }

  }
  plot(0,0,xlim=range(0,1),ylim=range(0,1),type="n",yaxs="i",xaxt="n",yaxt="n",bty="n")
  if(!usemid) {
    points(0.25,0.25,pch=1)
    text(0.3,0.25,colnames(rats)[3],pos=4,cex=cex)
    points(0.25,0.5,pch=2)
    text(0.3,0.5,colnames(rats)[1],pos=4,cex=cex)
  }
  else {
    points(0.25,0.25,pch=1)
    text(0.3,0.25,colnames(rats)[4],pos=4,cex=cex)
    points(0.25,0.5,pch=2)
    text(0.3,0.5,colnames(rats)[2],pos=4,cex=cex)
  }
  
  ow <- options("warn")$warn
  options(warn=-1)
  par(old.par)
  options(warn=ow)
}

setGeneric("plot")

setMethod("plot",c("QCStats","missing"),function(x,y,...) plot.qc.stats(x,...))



#old, DEPRECATED, WILL GO SOON!

getTao <- function(name) {
  .Deprecated("qc.get.tao","simpleaffy","This function (getTao) is DEPRECATED and WILL disappear from future versions of simpleaffy. Please see help(\"simpleaffy-deprecated\") for more details.\n")
  setQCEnvironment(name);
  qc.get.tau()
}

getAlpha1 <- function(name) {
  .Deprecated("qc.get.alpha1","simpleaffy","This function (getAlpha1) is DEPRECATED and WILL disappear from future versions of simpleaffy. Please see help(\"simpleaffy-deprecated\") for more details.\n")
  setQCEnvironment(name);
  qc.get.alpha1()
}

getAlpha2 <- function(name) {
  .Deprecated("qc.get.alpha2","simpleaffy","This function (getAlpha2) is DEPRECATED and WILL disappear from future versions of simpleaffy. Please see help(\"simpleaffy-deprecated\") for more details.\n")
  setQCEnvironment(name);
  qc.get.alpha2()
}

getActin3 <- function(name) {
  .Deprecated("qc.get.ratios","simpleaffy","This function (getActin3) is DEPRECATED and WILL disappear from future versions of simpleaffy. Please see help(\"simpleaffy-deprecated\") for more details.\n")
  setQCEnvironment(name);
  rats <- qc.get.ratios()
  toks <- unlist(strsplit(names(rats),"/"))
  rats <- unlist(rats)
  j <- cbind(toks,rats)
  r <- (j[j[,1] == "actin3",2])[1]
  names(r) <- NULL
  r
}

getActinM <- function(name) {
  .Deprecated("qc.get.ratios","simpleaffy","This function (getActinM) is DEPRECATED and WILL disappear from future versions of simpleaffy. Please see help(\"simpleaffy-deprecated\") for more details.\n")
  setQCEnvironment(name);
  rats <- qc.get.ratios()
  toks <- unlist(strsplit(names(rats),"/"))
  rats <- unlist(rats)
  j <- cbind(toks,rats)
  r <- (j[j[,1] == "actinM",2])[1]
  names(r) <- NULL
  r
}

getActin5 <- function(name) {
  .Deprecated("qc.get.ratios","simpleaffy","This function (getActin5) is DEPRECATED and WILL disappear from future versions of simpleaffy. Please see help(\"simpleaffy-deprecated\") for more details.\n")
  setQCEnvironment(name);
  rats <- qc.get.ratios()
  toks <- unlist(strsplit(names(rats),"/"))
  rats <- unlist(rats)
  j <- cbind(toks,rats)
  r <- (j[j[,1] == "actin5",2])[1]
  names(r) <- NULL
  r
}

getGapdh3 <- function(name) {
  .Deprecated("qc.get.ratios","simpleaffy","This function (getGapdh3) is DEPRECATED and WILL disappear from future versions of simpleaffy. Please see help(\"simpleaffy-deprecated\") for more details.\n")
  setQCEnvironment(name);
  rats <- qc.get.ratios()
  toks <- unlist(strsplit(names(rats),"/"))
  rats <- unlist(rats)
  j <- cbind(toks,rats)
  r <- (j[j[,1] == "gapdh3",2])[1]
  names(r) <- NULL
  r
}

getGapdhM <- function(name) {
  .Deprecated("qc.get.ratios","simpleaffy","This function (getGapdhM) is DEPRECATED and WILL disappear from future versions of simpleaffy. Please see help(\"simpleaffy-deprecated\") for more details.\n")
  setQCEnvironment(name);
  rats <- qc.get.ratios()
  toks <- unlist(strsplit(names(rats),"/"))
  rats <- unlist(rats)
  j <- cbind(toks,rats)
  r <- (j[j[,1] == "gapdhM",2])[1]
  names(r) <- NULL
  r
}

getGapdh5 <- function(name) {
  .Deprecated("qc.get.ratios","simpleaffy","This function (getGapdh5) is DEPRECATED and WILL disappear from future versions of simpleaffy. Please see help(\"simpleaffy-deprecated\") for more details.\n")
  setQCEnvironment(name);
  rats <- qc.get.ratios()
  toks <- unlist(strsplit(names(rats),"/"))
  rats <- unlist(rats)
  j <- cbind(toks,rats)
  r <- (j[j[,1] == "gapdh5",2])[1]
  names(r) <- NULL
  r
}

getAllQCProbes <- function(name) {
  .Deprecated("qc.get.probes","simpleaffy","This function (getAllQCProbes) is DEPRECATED and WILL disappear from future versions of simpleaffy. Please see help(\"simpleaffy-deprecated\") for more details.\n")
  setQCEnvironment(name)
  qc.get.probes()
}

getBioB <- function(name) {
  .Deprecated("qc.get.spikes","simpleaffy","This function (getBioB) is DEPRECATED and WILL disappear from future versions of simpleaffy. Please see help(\"simpleaffy-deprecated\") for more details.\n")
  setQCEnvironment(name)
  qc.get.spikes()["bioB"]
}


getBioC <- function(name) {
  .Deprecated("qc.get.spikes","simpleaffy","This function (getBioC) is DEPRECATED and WILL disappear from future versions of simpleaffy. Please see help(\"simpleaffy-deprecated\") for more details.\n")
  setQCEnvironment(name)
  qc.get.spikes()["bioC"]
}


getBioD <- function(name) {
  .Deprecated("qc.get.spikes","simpleaffy","This function (getBioD) is DEPRECATED and WILL disappear from future versions of simpleaffy. Please see help(\"simpleaffy-deprecated\") for more details.\n")
  setQCEnvironment(name)
  qc.get.spikes()["bioD"]
}


getCreX <- function(name) {
  .Deprecated("qc.get.spikes","simpleaffy","This function (getCreX) is DEPRECATED and WILL disappear from future versions of simpleaffy. Please see help(\"simpleaffy-deprecated\") for more details.\n")
  setQCEnvironment(name)
  qc.get.spikes()["creX"]
}


getAllSpikeProbes <- function(name) {
  .Deprecated("qc.get.spikes","simpleaffy","This function (getAllSpikeProbes) is DEPRECATED and WILL disappear from future versions of simpleaffy. Please see help(\"simpleaffy-deprecated\") for more details.\n")
  setQCEnvironment(name)
  qc.get.spikes()
}


haveQCParams <- function(name) {
  .Deprecated("qc.get.spikes","simpleaffy","This function (haveQCParams) is DEPRECATED and WILL disappear from future versions of simpleaffy. Please see help(\"simpleaffy-deprecated\") for more details.\n")
  qc.have.params(name)
}


