.if.probeset.not.found <- function(x) {
  return("NoAnno")
}

.strip.list <- function(x) {
  mapply(function(y,z) {
                        if(length(y) > 1) {
                            warning(paste("'",z,"' has more than one entry in annotation list. Taking the first one."));
    	 		    print(y)
                        }
                        return(y[1]);
                      },x,names(x));
}


get.annotation <- function (x, cdfname,verbose=FALSE) {
    library(cdfname, character.only = TRUE)
    symb <- .strip.list(mget(x, envir = get(paste(cdfname, "SYMBOL",
        sep = "")),ifnotfound=list(.if.probeset.not.found)))
    desc <- .strip.list(mget(x, envir = get(paste(cdfname, "GENENAME",
        sep = "")),ifnotfound=list(.if.probeset.not.found)))
    accno <- .strip.list(mget(x, envir = get(paste(cdfname, "ACCNUM",
        sep = "")),ifnotfound=list(.if.probeset.not.found)))
    uni <- .strip.list(mget(x, envir = get(paste(cdfname, "UNIGENE",
        sep = "")),ifnotfound=list(.if.probeset.not.found)))
    ok <- (symb != "NoAnno") & (desc != "NoAnno") & (accno!= "NoAnno") & (uni != "NoAnno")
    names(ok) <- x
    if(!ok && verbose) { warning(paste("value for '",names(ok)[ok], "' not found", sep = ""), call. = FALSE) }
    acc.lnk <- paste("=HYPERLINK(\"http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=search&db=nucleotide&term=",
        accno, "\",\"", accno, "\")", sep = "")
    acc.lnk[!ok] <- "NoAnno"
    uni.lnk <- paste("=HYPERLINK(\"http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=search&db=unigene&term=",
        uni, "&dopt=unigene\",\"", uni, "\")", sep = "")
    uni.lnk[!ok] <- "NoAnno"
    res <- cbind(symb, acc.lnk, uni.lnk, desc)
    res[res=="NoAnno"] <- "No Annotation Found"
    colnames(res) <- c("gene name", "accession", "unigene", "description")
    return(res)
}


write.annotation <- function(summary,file="results/annotation.table.xls") {
  write.table(summary,file=file,sep="\t",quote=F,col.names=NA)
}

results.summary <- function(results,cdfname) {
  res <- cbind(means(results),fc(results),sapply(fc(results),function(x) { if(x <0) { -1 * 2 ^ (-1 * x) } else { 2^x } }),tt(results),get.annotation(names(fc(results)),cdfname));
  cns <- colnames(res);
  cns[3]<-"log2(fc)";
  cns[4]<-"fc";
  cns[5]<-"tt";
  colnames(res) <- cns;
  return(res);
}


journalpng <- function(file="figure.png",width=4, height=4,res=300) {
  bitmap(file=file,type="png16m",width=width,height=height,res=res)
}

screenpng <- function(file="figure.png",width=4, height=4,res=72) {
  bitmap(file=file,type="png16m",width=width,height=height,res=res)
}
