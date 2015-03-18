\name{readDepth-package}
\docType{package}
\title{
  What the package does (short line)
}
\description{
  Utilizes multi-core systems with the help of the foreach and doMC
  packages, degrades gracefully to single core systems.
}
\details{
  \tabular{ll}{
    Package: \tab readDepth\cr
    Type: \tab Package\cr
    Version: \tab 1.0\cr
    Date: \tab 2010-05-12\cr
    License: \tab MIT\cr
    LazyLoad: \tab yes\cr
  }
  

  ~~ An overview of how to use the package, including the most important ~~
  ~~ functions ~~
}
\author{
  Chris Miller
  
  Maintainer: Chris Miller <chrisamiller@gmail.com>
}
\references{
  TBA
}
\keyword{ package }
\examples{

  ##create a rdObject, set bin params
  # x <- new("rdObject")

  ##fill the bins with reads
  # x.reads <- readDepth(x)

  ##do mapability correction
  # x.mapCorr <- rd.mapCorr(x.reads)

  ##do gc correction on the read counts
  # x.gcCorr <- rd.gcCorrect(x.mapCorr)

  ##return regions of CN gain and loss
  # segs <- rd.cnSegment(x.gcCorr)

}
