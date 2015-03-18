##-----------------------------------------------
## takes a breakpoint file in bed format and matches CNA segment 
## ends with the breakpoints if they're within binSize/2. It then
## adjusts segment edges if there is exactly one breakpoint 
## that matches.

rd.matchBreakpoints <- function(rdo,segs){

  ## read in bps from file
  bps <- read.table(paste("annotations/breakpoints.dat",sep=""), sep="\t")

  ##add an index  
  segs = cbind(segs,data.frame(index=1:length(segs[,1])))

  ##get only the alterations
  altSegs = subset(segs,(seg.mean > rdo@binParams$gainThresh/(rdo@binParams$med/2) | seg.mean < rdo@binParams$lossThresh/(rdo@binParams$med/2)))

  #count up total length for stderr ouput 
  bplen <<- c()
  
  ##match them up with the edges of segments
  corrSegs = foreach(chr=names(which(table(altSegs$chrom)!=0)), .combine="rbind") %do% {
    breakPointInt(rdo, bps[which(bps$V1==chr),], altSegs[which(altSegs$chrom == chr),], chr) 
  }

  ## we get arithmetic/rounding errors on
  ## some architectures unless we do this first
  corrSegs$loc.start = round(corrSegs$loc.start)
  corrSegs$loc.end = round(corrSegs$loc.end)
  segs$loc.start = round(segs$loc.start)
  segs$loc.end = round(segs$loc.end)

  ## put the altered ones back into the full data frame
  for(i in 1:length(corrSegs[,1])){
    ##get position of row to replace
    #which(row.names(segs)==row.names(corrSegs[i,]))
    pos = which(segs$index == corrSegs[i,]$index)    
    segs[pos,]=corrSegs[i,]

    ##adjust windows ahead and behind to compensate for change
    if(pos!=1){
      if(segs[pos-1,]$chrom == segs[pos,]$chrom){
        if(segs[pos-1,]$loc.end+1 != segs[pos,]$loc.start){
          segs[pos-1,]$loc.end = segs[pos,]$loc.start-1  
        }
      }
    }
    if(pos!=max(segs$index)){
      if(segs[pos+1,]$chrom == segs[pos,]$chrom){
        if(segs[pos+1,]$loc.start-1 != segs[pos,]$loc.end){
          segs[pos+1,]$loc.start = segs[pos,]$loc.end+1
        }     
      }
    }
  }

  if(verbose){
    cat("total number of adjusted breakpoints:",length(bplen),"\n")
    cat("average resolution of adjusted breakpoints:",mean(bplen),"\n")
  }

  return(segs[,1:5]) 
}

##---------------------------------------------
## actually do the intersection
breakPointInt <- function(rdo, bps, subSeg, chr){

  ##get a list of bps
  bpList = c()
  pos = c()
  ##exclude beginning and end of chr - not bps
  z = sort(append(subSeg$loc.start[1:length(subSeg$loc.start)],subSeg$loc.end[1:(length(subSeg$loc.start))]))
  for(i in 1:(length(z)-1)){
    if(!(z[i]+1 == z[i+1])){
      if(!(z[i] == 1) && !(z[i] == getChrLength(chr,rdo@entrypoints))){
        bpList = append(bpList,z[i])
        pos = append(pos,i)
      }
    }
  }
  bpList = append(bpList,z[length(z)])
  
  ##create a window around each bp
  segRange <- IRanges(start <- bpList-rdo@binParams$binSize/2, end <- bpList+rdo@binParams$binSize/2)  
  bpRange <- IRanges(start <- bps[,2], end <- bps[,3])
  
  reg <- as.matrix(findOverlaps(segRange,bpRange))
  singles = as.numeric(names(which(table(reg[,1])==1)))

  if(verbose){
    cat(chr,":\t",length(singles),"breakpoints being adjusted \n")
  }

  ##do adjustments to breakpoints based on paired-end data
  if(!(length(singles)==0)){
    for(i in 1:length(singles)){
      ##get breakpoint position
      bp = bps[reg[which(reg[,1]==singles[i]),2],]
      bplen <<- append(bplen,(bp$V3-bp$V2))


      ##find proper segment
      segPos = which(subSeg$loc.start==start(segRange)[reg[which(reg[,1]==singles[i])]]+rdo@binParams$binSize/2)
      if(!(length(segPos)==0)){
        subSeg[segPos,]$loc.start = (bp$V2+bp$V3)/2
      }
      segPos = which(subSeg$loc.end==start(segRange)[reg[which(reg[,1]==singles[i])]]+rdo@binParams$binSize/2)
      if(!(length(segPos)==0)){
        subSeg[segPos,]$loc.end = (bp$V2+bp$V3)/2
      }
    }
  }
  return(subSeg)
}
