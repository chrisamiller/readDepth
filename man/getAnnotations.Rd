\name{getAnnotations}
\alias{getAnnotations}
\title{
  getAnnotations
}
\description{
  connect to the google code repository and retrieve the 
  appropriate annotation files
}
\usage{
  getAnnotations(readLength, sex, genome="hg18", bs=FALSE)
}
\arguments{
  \item{readLength}{
    the average length of the reads, in bp.  
  }
  \item{sex}{
    either "male" or "female"
  }
  \item{genome}{
    genome build, in the form of "hg18"
  }
  \item{bs}{
    whether the reads are the result of bisulfite sequencing
    (TRUE/FALSE)
  }
}
\details{
  As of v1.0, the repo contains info for both hg18 and hg19 
  with read lengths of: 27, 36, 55, 76, 100, and 200.
 
  Bisulfite annotations have only been calculated 76bp reads on hg18

  Instructions for creating your own annotation files for parameter
  sets not listed here are available at: 
  http://code.google.com/p/readdepth
}
\examples{
  #  library(readDepth)
  #  getAnnotations(readLenght=76, sex="male", genome="hg19")
}
