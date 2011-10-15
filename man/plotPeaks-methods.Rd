\name{plotPeaks}
\docType{methods}
\alias{plotPeaks}
\alias{plotPeaks-methods}
\alias{plotPeaks,CNAnorm-method}
\title{ Methods for Function plotPeaks in Package `CNAnorm' }
\description{
\code{plotPeaks} plot annotated distribution of ratio Test/Normal 
}

\usage{
\S4method{plotPeaks}{CNAnorm}(object, special1 = character(0), special2 = character(0),
    show ='suggested', adjust = get.adjust(object), n = get.n(object), ...)
}

\arguments{
 
\item{object}{An object of Class \code{"CNAnorm"}}
\item{special1}{The chromosome/contig whose distribution will be shown with a
    different color }
\item{special2}{The chromosome/contig whose distribution will be shown with a
    different color }
\item{show}{A character verctor with one or both of the following:
    \code{"suggested"}, \code{"validated"}. Specify what need to be plotted
}
\item{adjust}{ The adjust parameter for function \code{"density"}. Default will 
retreive the values use for \code{peakPloidy}}
\item{n}{ The n parameter for function \code{"density"}. Default will 
retreive the values use for \code{peakPloidy}}
\item{...}{Further arguments to pass to the function \code{plot}}
        
}


\author{Stefano Berri \email{s.berri@leeds.ac.uk}}

\examples{
data(CN)
# see peakPloidy documentation to know how object CN is created
plotPeaks(CN, special1 = 'chrX', special2 = 'chrY')
}
\seealso{
\code{\link{plot}}, \code{\link{validation}}, \code{\link{peakPloidy}}
}

