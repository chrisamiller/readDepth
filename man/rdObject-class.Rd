\name{rdObject-class}
\Rdversion{1.0}
\docType{class}
\alias{rdObject-class}

\title{Class "rdObject" ~~~ }

\description{
Objects of this class hold data used by the readDepth package, 
including input parameters, depth-of-coverage values, entrypoint information,
etc.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("rdObject")}.

In order to successfully create this class, you'll need to have the
proper directory structure set up. See the readDepth documentation 
for details (http://code.google.com/p/readdepth)

}

\section{Slots}{
  \describe{
    \item{\code{params}:}{Object of class \code{"numeric"}
       a named list containing input parameters
    }
    \item{\code{binParams}:}{Object of class \code{"data.frame"}
       a data frame that contains calculated values, like thresholds, 
       mean # of reads per bin, ploidy calculations, etc.
    }
    \item{\code{entrypoints}:}{Object of class \code{"data.frame"}
       a data frame containing the entrypoints and mapability values
    }
    \item{\code{readInfo}:}{Object of class \code{"data.frame"}
       a data frame containing information on the input data (reads)
    }
    \item{\code{chrs}:}{Object of class \code{"list"} 
       a named list with each slot containing a dataframe that holds
       a given chromosome's read depth, gc content, and mapability 
       for each bin.
    }
  }
}
\section{Methods}{
No methods defined with class "rdObject" in the signature.
}
\references{ readDepth: a parallel R package for detecting copy-number alterations from short sequencing reads (in preparation)}
\author{ Chris Miller (chrisamiller@gmail.com) }
\keyword{classes}
