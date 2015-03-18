\name{trimSegmentEnds}
\alias{trimSegmentEnds}
\title{
  trimSegmentEnds
}
\description{
  trim the ends of segments such that they're called
  normal (CN = 2.0) until the first bin that isn't NA
  
  useful for chromosomes like 22, where there are tens of
  megabases that are unmappable at beginning. Without this,
  those regions would be called with the same CN as the first
  few mappable bins.
  
}
\usage{
  trimSegmentEnds(segs, rdo)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{segs}{
    a data frame that is the result of calling rd.cnSegments
  }
  \item{rdo}{
    a readDepth object created with new("rdObject")
  }
}
\value{
  the segment data frame adjusted to have the ends trimmed
}
\examples{

  ##get regions of CN gain and loss
  #  segs <- rd.cnSegments(x)

  ##trim the segment ends 
  #  trimSegmentEnds(segs, rdo)
}
