##------------------------------------------------
## Given an rd object, add the number of reads in
## each bin. 
readDepth <- function(rdo){

  ##set multicore options for reading if specified
  if("readCores" %in% names(rdo@params)){
    options(cores = as.numeric(rdo@params[["readCores"]]))
  }

  if(verbose){
    cat("Binning Start:",date(),"\n")
  }

  ##do the binning for each chr in parallel
  mcoptions <- list(preschedule = FALSE)
  b = foreach(chr=rdo@entrypoints$chr, .combine="combineBins",.options.multicore=mcoptions) %dopar% {
    doBinning(rdo@params, rdo@binParams, rdo@entrypoints, rdo@readInfo, chr)
  }
  rdo@chrs = b
  closeAllConnections()
  gc() #just in case


  if(verbose){
    cat("Binning End:",date(),"\n")
    cat("-------------------------------\n")
    #plot distribution
    plotDist(rdo, filename="output/dist.rawreads.pdf")
  }

  ##reset multicore options
  if("maxCores" %in% names(rdo@params)){
    options(cores = as.numeric(rdo@params[["maxCores"]]))
  }

  return(rdo)
}



##--------------------------------------------------
## get the read depth for this chromosome and
## output a data frame in the proper position of the list
##
doBinning <- function(params, binParams, entrypoints, readInfo, chr){
  pref = paste(chr,":\t",sep="")

  if(verbose){
    cat(pref,"Starting\n")
  }
  
  binSize=binParams$binSize
  chrLen = getChrLength(chr,entrypoints)

  ## chunk the file if necessary to avoid running out of memory
  filename <- paste(chr,".bed",sep="")   
  filepath <- paste("reads/",chr,".bed",sep="")  	 
  fileLen <- readInfo[which(readInfo$file==filepath),]$len  

  b <- foreach(pos=seq(0,fileLen,by=params[["chunkSize"]]), .combine="+") %do% {
    binReads(binSize=binSize, chr=chr, chrLen=chrLen, filename=filepath, start=pos, chunkSize=params[["chunkSize"]], zfileLen=fileLen)
  }
  if(verbose){
    cat(pref,"Done \n")
  }

  return(listify(data.frame(rd=b),chr))
}


##--------------------------------------------------
## map reads to bins for a given file
##
binReads <- function(binSize,chr,chrLen,filename,start,chunkSize,zfileLen){
  ## create an IRange for bins
  bins <- IRanges(start = (0:(ceiling(chrLen/binSize)-1)*binSize)+1, end = (1:ceiling(chrLen/binSize))*binSize)
  end(bins[length(bins)]) <- chrLen
  cat(chr,": reading chunk",(round(start/chunkSize)+1),"of",ceiling(zfileLen/chunkSize),"\n")
  
  ## input the reads
  rawreads = scan(pipe(paste("cut -f2",filename)), what=0, sep="\t", skip=start, nlines=chunkSize, quiet=TRUE)
  reads = IRanges(start = rawreads, end=rawreads)
  ##free up some memory
  fullreads = NULL
  ##figure out which reads fall in which bins
  binnedReads <- as.matrix(findOverlaps(bins,reads))

  nonZero = table(binnedReads[,1])

  a=rep(0,length(bins))
  b=names(nonZero)
  for(i in 1:length(b)){
    a[as.numeric(b[i])]=nonZero[i]
  }
  return(a)
}

