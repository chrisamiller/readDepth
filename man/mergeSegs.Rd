\name{mergeSegs}
\alias{mergeSegs}
\title{
  mergeSegs
}
\description{
   merges segments if they are adjacent and both gains or losses
   (useful for some types of plotting)
}
\usage{
  mergeSegs(segs,rdo)
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
  a data frame containing merged segments of CN gain and loss and all
  intervening normal segments
}
\seealso{
  'rd.cnSegments'
}

\examples{

  ##get regions of CN gain and loss
  #  segs <- rd.cnSegments(x)

  ##extract just the gains and losses
  #  mergedSegs <- mergeSegs(segs,x)
}
