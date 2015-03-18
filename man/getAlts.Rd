\name{getAlts}
\alias{getAlts}
\title{
  getAlts
}
\description{
  takes a data frame of segments and extracts just those containing
  copy-number gains and losses (that exceed the gain/loss thresholds)
}
\usage{
  getAlts(segs,rdo)
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
  a data frame containing only segments of CN gain and loss
}
\seealso{
  'rd.cnSegments'
}

\examples{

  ##get regions of CN gain and loss
  #  segs <- rd.cnSegments(x)

  ##extract just the gains and losses
  #  alts <- getAlts(segs,x)
}
