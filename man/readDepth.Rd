\name{readDepth}
\alias{readDepth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  readDepth
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  Counts the number of sequence reads that fall into identically
  sized bins across the entire genome. 
}
\usage{
  readDepth(rdo)
}
\arguments{
  \item{rdo}{
    \code{rdo}: a readDepth Object, created with the command 'new("rdObject")'
  }
}
\details{  
  %%  ~~ If necessary, more details than the description above ~~
  Assumes that the correct directory structure is in place (see
  full description in readDepth documentation)

  In brief, the directory "reads" should be present and should contain
  read positions in bed format.
}
\value{
  
  returns the original rdObject with the @chrs attribute updated to
  contain a list of N data frames (one for each chromosome). Each data
  frame contains one column "rd" containing the counts for each bin.
}

\examples{

  ##create a rdObject, set bin params
  # x <- new("rdObject")

  ##fill the bins with reads
  # x.reads <- readDepth(x)

}
