##-------------------------------------------------
## Common functions used by the readDepth package
##


##-------------------------------------------------
## load libraries

library('methods')
library('foreach')
library('IRanges')
library('doMC')


##-------------------------------------------------
## read in entrypoints file, return a data frame
##
readEntrypoints <- function(){
  p=read.table("annotations/entrypoints",sep="\t",quote="",colClasses=c("character","numeric","numeric"))
  names(p)=c("chr","length","ploidy")  
  return(p)  
}

##------------------------------------------------
##strip leading and trailing whitespace
##
trim <- function (x) gsub("^\\s+|\\s+$", "", x)


##--------------------------------------------------
##  Do a wc on each read file to find out how many reads there 
##  are. This is way faster than counting with R. Eventually,
##  we can write a little C function to do this so we don't  
##  require a shell with wc installed (allow windows port)
##
getReadInfo <- function(){
  ## check to see if they've already been counted
  lenFileName="reads/readInfo"
  len = 0
  
  if(file.exists(lenFileName)){
    flist = read.table(lenFileName)
    names(flist) = c("len","file")
  } else {
    lenList = foreach(e=Sys.glob("reads/*.bed"), .combine="append") %do%{
      strsplit(trim(system(paste('wc -l ', e, sep=""), intern=T)),"\\s+",perl=T)[[1]]
    }
    flist = data.frame(len=as.numeric(lenList[seq(1,length(lenList)-1,2)]), file=lenList[seq(2,length(lenList),2)])
    write.table(flist,file=lenFileName,sep="\t",quote=FALSE, row.names=FALSE, col.names=FALSE)
  }  
  return(flist)  
}


##--------------------------------------------------
##  make sure each entrypoint has a corresponding bed file
##
verifyFiles <- function(chrs){
  for(i in 1:length(chrs)){
    filename = paste("reads/",chrs[i],".bed",sep="")
    if(!(file.exists(filename))){
      stop("file: '",filename,"' does not exist\n")
    }
  }
}


##-------------------------------------------------
## sum up the already-stored # of reads
##
getNumReads <- function(flist){
  len = sum(as.numeric(flist$len))
  return(len)
}

##--------------------------------------------------
## generate a distribution by modelling the overdispersed
## poisson with the negative binomial distribution
## d = var/mean
##
rpois.od<-function (n, lambda,d=1) {
  if (d==1)
    rpois(n, lambda)
  else
    rnbinom(n, size=(lambda/(d-1)), mu=lambda)
}


##----------------------------------------------------
## return the length or ploidy for a given chr identifier
## 
getChrLength <- function(chr,entries){
  return(entries[which(entries$chr==chr),]$length)
}

getChrPloidy <- function(chr,entries){
  return(entries[which(entries$chr==chr),]$ploidy)
}


##----------------------------------------------------
## create a list of the appropriate size and place
## input in the appropriate bin 
##
listify <- function(input,chr){
  lst = vector("list")
  lst[[chr]] = input
  return(lst)
}


##----------------------------------------------------
## function to combine lists returned from parallel
## intersection procedure
##
combineBins <- function(a,b){ 
  for(i in 1:length(a)){
    name=names(a)[i]
    if(is.null(b$name)){
      b[[name]] = a[[i]]
    }
  }
  return(b)
}
 

##---------------------------------------------------
## convert chromosome number to chr name
##
chrName <- function(num){
  if(num == 23){ num <- "X"}
  if(num == 24){ num <- "Y"}
  
  return(paste("chr",num,sep=""))  
}


##---------------------------------------------------
## sum the total lengths of the annotations in a 
## bed file. This is equal to coverage if annotations
## are non-overlapping
##
bedAnnotationLength <- function(e){
  if(file.exists(e)){
    f=gzfile(e)
    a=scan(f,what=0,quiet=TRUE)
    close(f)
    return(sum((a[seq(2,(length(a)),2)]-a[seq(1,(length(a)-1),2)]+1)))
  }
  return(0)
}


##--------------------------------------------------
## convert read depth to log2 value,
## based on median
##
logScore <- function(val,med){
  if(is.na(val)){
    return(NA)

  }else if(val<=0){
    return(floor(log2(1/med)))

  }else{
    return(log2(val/med))
  }    
}


##----------------------------------------------
## given a set of parameters, plot the peaks
## and thresholds
##
plotWindows <- function(windowSize, genomeSize, divGain, divLoss, fdr, numReads, oDisp, ploidyPerc, med){
  numWinds <- genomeSize/windowSize
  
  ##generate distributions
  d2 <- rpois.od(numWinds*ploidyPerc$diploidPerc,med,oDisp)
  d3 <- rpois.od(numWinds*ploidyPerc$triploidPerc,(med*1.5),oDisp)
  d1 <- rpois.od(numWinds*ploidyPerc$haploidPerc,(med*0.5),oDisp)

  hrange=10*sqrt(med*oDisp)
  hist(d2,breaks=seq(-100000,100000,20),xlim=c(med-hrange,med+hrange),main=paste("Window Size: ",windowSize,sep=""))
  mtext(paste("FDR: ",round(fdr,5),sep=""))
  hist(d3,breaks=seq(-100000,100000,20),add=T,col="red")
  hist(d1,breaks=seq(-100000,100000,20),add=T,col="red")
  abline(v=med,col="blue")
  abline(v=med*(0.5),col="blue")
  abline(v=med*(1.5),col="blue")
  abline(v=divGain,col="green")
  abline(v=divLoss,col="green")
}

##----------------------------------------------
## plot the segments for a given chromosome
##
plotSegs <- function(rdo,segs,chr){

  st = 1
  sp = rdo@entrypoints[which(rdo@entrypoints$chr == chr),]$length
  
  winds = rdo@chrs[[chr]]$rd
  #print(winds)
  binSize = rdo@binParams$binSize  
  pos = seq(binSize/2,(((length(winds)-1)*binSize)-binSize/2),binSize)
  pos = append(pos,sp)
  
  par(mar=c(5, 4, 4, 4) + 0.1)
  plot(pos,winds,ylab="number of reads", xlab="position (bp)")

  abline(h=rdo@binParams$med,col="blue")
  abline(h=rdo@binParams$gainThresh,col="green")
  abline(h=rdo@binParams$lossThresh,col="green")

  asegs = segs[which(segs$chrom == chr),]
  for(i in 1:length(asegs[,1])){
    lines(c(asegs[i,2],asegs[i,3]), c(asegs[i,5],asegs[i,5]), col="red",lwd=3)
  }
  
  par(new=T)
  plot(-10000,-10000,ylim=c(0,(max(winds,na.rm=TRUE)/rdo@binParams$med)), xlim=c(1,sp),axes=F,xlab="", ylab="")
  axis(4, ylim=c(0,max(winds,na.rm=TRUE)/rdo@binParams$med), col="red",col.axis="red")
  mtext("Copy Number",side=4,col="red",line=2.5)

}


##--------------------------------------------------
## plot overlapping histograms
##
plotOverlappingHist <- function(a, b, colors=c("white","gray20","gray50"),
                                breaks=NULL, xlim=NULL, ylim=NULL){
  
  ahist=NULL
  bhist=NULL

  if(!(is.null(breaks))){
    ahist=hist(a,breaks=breaks,plot=F)
    bhist=hist(b,breaks=breaks,plot=F)
  } else {
    ahist=hist(a,plot=F)
    bhist=hist(b,plot=F)

    dist = ahist$breaks[2]-ahist$breaks[1]
    breaks = seq(min(ahist$breaks,bhist$breaks),max(ahist$breaks,bhist$breaks),dist)

    ahist=hist(a,breaks=breaks,plot=F)
    bhist=hist(b,breaks=breaks,plot=F)
  }

  if(is.null(xlim)){
    xlim = c(min(ahist$breaks,bhist$breaks),max(ahist$breaks,bhist$breaks))
  }

  if(is.null(ylim)){
    ylim = c(0,max(ahist$counts,bhist$counts))
  }

  overlap = ahist
  for(i in 1:length(overlap$counts)){
    if(ahist$counts[i] > 0 & bhist$counts[i] > 0){
      overlap$counts[i] = min(ahist$counts[i],bhist$counts[i])
    } else {
      overlap$counts[i] = 0
    }
  }

  plot(ahist, xlim=xlim, ylim=ylim, col=colors[1])
  plot(bhist, xlim=xlim, ylim=ylim, col=colors[2], add=T)
  plot(overlap, xlim=xlim, ylim=ylim, col=colors[3], add=T)
}


##--------------------------------------------------
## plot actual and expected histograms
##
plotDist <- function(rdo,xmax=NULL,filename="output/hist.pdf",windSize=NULL){

  bins=c()
  for(i in 1:length(names(rdo@chrs))){
    bins = append(bins,rdo@chrs[[names(rdo@chrs)[i]]]$rd)
  }

  if(is.null(xmax)){
    xmax=rdo@binParams$gainThresh*2 ##max(bins,na.rm=TRUE) ##sort(bins[round(length(bins)*0.999)])
  }
  
  if(is.null(windSize)){
    windSize = round(xmax/100)
  }

  len=length(bins)
  e=rdo@entrypoints
  breaks=seq(0,xmax+windSize,windSize)
  oDisp=rdo@params["overDispersion"]

  med=rdo@binParams$med 

  p1 = rpois.od(len*rdo@binParams$hapPerc,med/2,oDisp)
  p2=rpois.od(len*rdo@binParams$dipPerc,med,oDisp)
  p3 = rpois.od(len*rdo@binParams$tripPerc,(med*(3/2)),oDisp)
  
  p=append(append(p1,p2),p3)

  bins = bins[which(bins<xmax & bins>0)]
  p = p[which(p<xmax & p>0)]

  pdf(file=filename)
  plotOverlappingHist(a=bins,b=p,breaks=breaks,xlim=c(0,xmax))
  abline(v=rdo@binParams$gainThresh,col="red")
  abline(v=rdo@binParams$lossThresh,col="red")
  dev.off()
}


##--------------------------------------------------
## plot three actual and expected histograms 
## usually used with raw read depth, post mapability
## correction, and post gc content correction
##
tripDist <- function(xmax=NULL,filename="output/hist.pdf",windSize=20, rdo2, rdo3, rdo4){
  pdf(file=filename, width=12,height=4)
  par(mfcol=c(1,3))
  
  for(rdo in c(rdo2,rdo3,rdo4)){
    bins=c()
    for(i in 1:length(names(rdo@chrs))){
      bins = append(bins,rdo@chrs[[names(rdo@chrs)[i]]]$rd)
    }

    if(is.null(xmax)){
      xmax=rdo2@binParams$gainThresh*2 ##max(bins,na.rm=TRUE) ##bins[round(length(bins)*0.999)]
    }
    len=length(bins)
    e=rdo@entrypoints
    breaks=seq(0,xmax,windSize)
    oDisp=rdo@params["overDispersion"]
    
    med=rdo@binParams$med
    
    p1 = rpois.od(len*rdo@binParams$hapPerc,med/2,oDisp)
    p2=rpois.od(len*rdo@binParams$dipPerc,med,oDisp)
    p3 = rpois.od(len*rdo@binParams$tripPerc,(med*(3/2)),oDisp)
    
    p=append(append(p1,p2),p3)
    
    bins = bins[which(bins<xmax & bins>0)]
    p = p[which(p<xmax & p>0)]    

    plotOverlappingHist(a=bins,b=p,breaks=breaks)
    abline(v=rdo@binParams$gainThresh,col="red")
    abline(v=rdo@binParams$lossThresh,col="red")
  }

  dev.off()
}


##-----------------------------------------------------------
## some simple output functions
##
writeThresholds <- function(rdo){
  a1=data.frame(gainThresh=(rdo@binParams$gainThresh/rdo@binParams$med)*2)
  a2=data.frame(lossThresh=(rdo@binParams$lossThresh/rdo@binParams$med)*2)
  a3=data.frame(binSize=rdo@binParams$binSize)

  write.table(t(cbind(a1,a2,a3)),file="output/thresholds.dat",sep="\t",quote=F,row.names=T,col.names=F)
}

writeSegs <- function(segs){
  write.table(segs, file="output/segs.dat", sep="\t", quote=F, row.names=F, col.names=F)
}

writeAlts <- function(segs,rdo){
  write.table(getAlts(segs,rdo), file="output/alts.dat", sep="\t", quote=F, row.names=F, col.names=F)
}


##------------------------------------------------------------
## Code for pulling down the annotation files 
## and sticking them in the right place

getAnnotations <- function(readLength, sex, genome="hg18", bs=FALSE){
  if(sex!="male" & sex!="female" & sex!="autosomes"){
    print("'sex' must be either \"male\", \"female\", or \"autosomes\"")
    return(0)
  }

  
  dlAnnFile <- function(url,outfile){
    download.file(url, outfile, quiet=TRUE)
    if(file.exists(outfile)){
      print(paste(url,"downloaded successfully"))
      return(1)
    } else {
      print(paste("Oops! Couldn't fetch ",url))
      print("Either you specified a read length or genome build that we don't have annotation")
      print("files prepared for, or your network connection isn't working.")
      print("Check the documentation at http://github.com/chrisamiller/readDepth/ for more information")
      return(0)
    }
  }

  entryurl=paste("https://xfer.genome.wustl.edu/gxfer1/project/cancer-genomics/readDepth/entrypoints.",genome,".",sex,sep="")
  entryfile=paste("annotations/entrypoints")    
  
  if(bs){ #bisulfite
    mapurl = paste("https://xfer.genome.wustl.edu/gxfer1/project/cancer-genomics/readDepth/mapability.bs.readLength",readLength,".",genome,".tar",sep="")
    mapfile = paste("annotations/mapability.bs.readLength",readLength,".",genome,".tar",sep="")
    gcurl = paste("https://xfer.genome.wustl.edu/gxfer1/project/cancer-genomics/readDepth/gcWinds.bs.readLength",readLength,".",genome,".tar",sep="")
    gcfile = paste("annotations/gcWinds.bs.readLength",readLength,".",genome,".tar",sep="")    
  } else { #normal
    mapurl = paste("https://xfer.genome.wustl.edu/gxfer1/project/cancer-genomics/readDepth/mapability.readLength",readLength,".",genome,".tar",sep="")
    mapfile = paste("annotations/mapability.readLength",readLength,".",genome,".tar",sep="")
    gcurl = paste("https://xfer.genome.wustl.edu/gxfer1/project/cancer-genomics/readDepth/gcWinds.readLength",readLength,".",genome,".tar",sep="")
    gcfile = paste("annotations/gcWinds.readLength",readLength,".",genome,".tar",sep="")    
  }

  #download the files, untar them
  if(dlAnnFile(mapurl, mapfile)){
    system(paste("tar -xf",mapfile,"-C annotations/"))
  }
  
  if(dlAnnFile(gcurl, gcfile)){
    system(paste("tar -xf",gcfile,"-C annotations/"))
  }

  dlAnnFile(entryurl, entryfile)
#    system(paste("tar -xf",entryfile,"-C annotations/"))
#  }
}

