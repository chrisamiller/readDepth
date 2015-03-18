\name{rd.cnSegments}
\alias{rd.cnSegments}
\title{
  Identifies regions of copy-number alteration
}
\description{  
  
  uses Circular Binary Segmentation (as implemented in the DNAcopy
  package) to segment the genome and returns regions that contain CN
  gains or losses
}
\usage{
  rd.cnSegments(rdo, onlyAlts=FALSE, minWidth=2, alpha=0.01, rmGaps=TRUE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rdo}{
    a readDepth object created with new("rdObject") and filled with read
    counts using the readDepth() command
  }
  \item{onlyAlts}{
    whether the function should return just those segments that 
    exceed the threshold values for gain and loss The Default is FALSE,
    which returns all segments.
  }
  \item{minWidth}{
    the minimum number of consecutive aberrant bins required before making
    a segment call.  Default is two 
  }
  \item{alpha}{
    the significance level for the CBS to accept change-points. Default 0.01
  }
  \item{rmGaps}{
    Controls whether adjacent segments have their boundaries adjusted to remove
    gaps. If so, the new segment boundary is set to halfway between the existing
    boundary calls.
  }
}
\value{
  returns a dataframe with 5 columns:
  chr  start  stop  #probes  copynumber

}
\seealso{
  'readDepth'
}
\examples{

  #create a rdObject, set bin params
  #  x <- new("rdObject")

  #fill the bins with reads
  #  x.reads <- readDepth(x)

  #do mapability correction on the read counts
  #  x.map <- rd.mapCorrect(x.reads)

  #do gc correction on the read counts
  #  x.gcCorr <- rd.gcCorrect(x.map)

  #return regions of CN gain and loss
  #  segs <- rd.cnSegments(x.gcCorr)

}
