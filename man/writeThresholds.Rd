\name{writeThresholds}
\alias{writeThresholds}
\title{
  writeThresholds
}
\description{
  write out a tab-delimited table of the calculated
  gain/loss thresholds to the file "output/thresholds.dat"
}
\usage{
  writeThresholds(rdo)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rdo}{
    a readDepth object created with new("rdObject")
  }
}
\examples{

  ##get regions of CN gain and loss
  #  x <- new("rdObject")

  ##write out the table
  #  writeThresholds(x)
}
