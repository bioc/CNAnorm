\name{plotGenome}
\docType{methods}
\alias{plotGenome}
\alias{plotGenome-methods}
\alias{plotGenome,CNAnorm-method}
\title{ Methods for Function plotGenome in Package `CNAnorm' }
\description{
\code{plotGenome} plot normalized ratio and optionally segmented and/or smoothed 
    normalized ratio values in Package `CNAnorm'. It also shows annotation. 
}

\usage{
\S4method{plotGenome}{CNAnorm}(object, maxRatio = 8, minRatio = -1,
    superimpose = character(0),  supLineColor = character(0),
    supLineCex = character(0), numHorLables = 10, colorful = FALSE, ...)
}

\arguments{
 
\item{object}{An object of Class \code{"CNAnorm"}}
\item{maxRatio}{The maximum ratio to be shown on the plot. Values or ratio
    greater than maxRatio will be displayed as red triangulars}
\item{minRatio}{The minimum ratio to be shown on the plot. Values or ratio
    smaller than minRatio will be displayed as red triangulars}
\item{superimpose}{A character verctor with one or both of the following:
    \code{"smooth"}, \code{"DNACopy"}}
\item{supLineColor}{A three element character vector with colors to be used for
    \code{first superimposed line}, \code{second superimposed line}, 
    \code{normalized ratio dots}. If none is provided, defaults are:
    c("black", "cyan", "grey60") }
\item{supLineCex}{A two element character vector with \code{cex} valeus to be
    used for width of \code{first superimposed line} and \code{second
    superimposed line}. If none is provided, defaults are: c(0.5, 0.5) }
\item{numHorLables}{. Number of maximum horizontal lables. The function will try
    to annotate numHorLables so that they are approximately equally spaced.}
\item{colorful}{A switch to decide if the background dots representing the 
    ratio of each window should be gray or colored according their value in 
    relation to the peak closest to the median}
\item{...}{Further arguments to pass to the function \code{plot}}
        
}

\value{An object of class \code{"CNAnorm"}}

\author{Stefano Berri \email{s.berri@leeds.ac.uk} and 
    Arief Gusnanto \email{a.gusnanto@leeds.ac.uk}}
\examples{
data(CN)
# see peakPloidy documentation to know how object CN is created
CN <- addDNACopy(CN)
CN <- validation(CN)
CN <- discreteNorm(CN)
plotGenome(CN, superimpose = 'DNACopy')
}
\seealso{
\code{\link{plot}}, \code{\link{par}}, \code{\link{peakPloidy}}
}

