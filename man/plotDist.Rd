\name{plotDist}
\alias{plotDist}
\title{
  plotDist
}
\description{
  plot overlapping histograms of the model distribution and the actual 
  distribution of reads
}
\usage{
  plotDist(rdo,xmax,filename,windSize)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rdo}{
    a readDepth object created with new("rdObject")
  }
  \item{xmax}{
    the scale of the xaxis goes from 0 until this value (default 2000)
  }
  \item{filename}{
    the output file (default "output/hist.pdf")
  }
  \item{windSize}{
    the window size for the histogram to use (default 20)
  }
}
\value{
  a histogram that helps determine how well the model fits the data
}
\examples{

  ##get regions of CN gain and loss
  #  rdo <- readDepth(x)

  ##does the model fit well
  #  plotDist(rdo)

  ##does it fit better after mapability correction?
  #  rdo2 = rd.mapCorrect(rdo)
  #  plotDist(rdo2)
  
  ##adjust the graphical params for a different view
  #  plotDist(rdo2, 1000, "output/postMapCorrect.pdf", 10)
}
