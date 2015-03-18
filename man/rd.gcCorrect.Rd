\name{rd.gcCorrect}
\alias{rd.gcCorrect}
\title{
  readDepth: GC-content correction
}
\description{
  takes a readDepth object containing read counts and uses loess
  correction to account for GC content bias. Updates the @chrs rd column
  and adds a "gc" column containing that bin's GC content percentage
}
\usage{
  rd.gcCorrect(rdo, meth=FALSE, outlierPercentage=0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rdo}{
    a readDepth object created with new("rdObject") and filled with read
    counts using the readDepth() command
  }
  \item{meth}{
    if this data consists of bisulfite reads and a methylation map is
    present in annotations/methMap, gc correction will be done while
    taking into account the bisulfite treatment that converts C -> U 
    and methylated bases that were protected from this conversion
  }
  \item{outlierPercentage}{
    Removes the highest and lowest N% of reads before loess correction
    This prevents overfitting near the extremes, which loess is 
    susceptible to. Default value is 0.01    
  }
}
\details{
  Assumes that files containing gc content and unmappable regions are in
  place (see full readDepth documentation for proper directory structure).
}
\value{
  returns a rd object with the read depth values altered to reflect
  gc content composition-induced bias.
}
\seealso{
  'readDepth'
}

\examples{

  ##create a rdObject, set bin params
  #  x <- new("rdObject")

  ##fill the bins with reads
  #  x.reads <- readDepth(x)

  ##do gc correction on the read counts
  #  x.map <- rd.mapCorrect(x.reads)

  ##do gc correction on the read counts
  #  x.gcCorr = rd.gcCorrect(x.map)

}
