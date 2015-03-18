\name{rd.matchBreakpoints}
\alias{rd.matchBreakpoints}
\title{
  refine edges of CNAs by matching them with nearby breakpoints
}
\description{
  takes a readDepth object containing read counts and parameters, along
  with the results of CN segmentation, and refines the edges of the CN
  regions by matching them with nearby breakpoints. 
}
\usage{
  rd.matchBreakpoints(rdo,segs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rdo}{
    a readDepth object created with new("rdObject") and filled with read
    counts using the readDepth() command
  }
  \item{segs}{
    a data frame containing positions of copy number segments, as produced
    by the rd.cnSegment function
  }
}
\details{
   Assumes a breakpoint list (in bed format) is present in the annotations
   directory (see manual)

   If the region around the window contains more than one breakpoint, it is
   unable to determine which bp is correct and makes no adjustment
}
\value{
  returns a modified copy of the segs data frame that was fed in.
}
\seealso{
  'readDepth'
}

\examples{

  ##create a rdObject, set bin params
  #x <- new("rdObject")

  ##fill the bins with reads
  #x.reads <- readDepth(x)

  ##do gc correction on the read counts
  #x.gcCorr = rd.gcCorrect(x.reads)

  ## map

  ## cn seg

  ##x.matchBreakpoints(x,segs)
}
