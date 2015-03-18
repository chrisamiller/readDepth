\name{plotSegs}
\alias{plotSegs}
\title{
  plotSegs
}
\description{
  quick and dirty way of visualizing copy-number segments on a single
  chromosome
}
\usage{
  plotSegs(rdo,segs,chr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rdo}{
    a readDepth object created with new("rdObject")
  }
  \item{segs}{
    a data frame that is the result of calling rd.cnSegments
  }
  \item{chr}{
    the name of the chromosome to plot
  }
}
\value{
  a plot with position on the X-axis and segments of CN on the Y axes
}
\seealso{
  'rd.cnSegments'
}
\examples{

  ##get regions of CN gain and loss
  #  segs <- rd.cnSegments(x)

  ##extract just the gains and losses
  #  plotSegs(x, segs, "chr22")
}
