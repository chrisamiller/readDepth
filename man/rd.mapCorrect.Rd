\name{rd.mapCorrect}
\alias{rd.mapCorrect}
\title{
  Correct read depth for unmappable regions
}
\description{
  takes a readDepth object containing read counts and updates the rd
  counts, account for repetitive and unmappable genomic regions. Also
  adds a "map" column containing that bin's percentage mapability
}
\usage{
  rd.mapCorrect(rdo, minMapability=0.25)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rdo}{
    a readDepth object created with new("rdObject") and filled with read
    counts using the readDepth() command
  }
  \item{minMapability}{
    the minimum amount of mapability required. If the mapability of a
    region is below this level, the read depth is set to NA. 
  }
}
\details{
  Assumes that files containing unmappable regions are in place (see
  readDepth documentation for proper directory structure). 
}
\value{
  returns a readDepth object with the read depth values altered to
  remove bias due to mappability
}
\seealso{
  'readDepth'
}

\examples{
  ##create a rdObject, set bin params
  #  x <- new("rdObject")

  ##fill the bins with reads
  #  x.reads <- readDepth(x)

  ##do mapping correction on the read counts
  #  x.mapCorr = rd.mapCorrect(x.reads)
}
