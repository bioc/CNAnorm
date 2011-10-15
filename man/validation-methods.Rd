\name{validation}
\docType{methods}
\alias{validation}
\alias{validation-methods}
\alias{validation,CNAnorm-method}
\title{ Methods for Function addSmooth in Package `CNAnorm' }
\description{
\code{validation} segment and smooth perform ratio values in Package `CNAnorm' 
}

\usage{
\S4method{validation}{CNAnorm}(object, peaks = sugg.peaks(object), 
    ploidy = sugg.ploidy(object))
}

\arguments{
 
\item{object}{An object of Class \code{"CNAnorm"}}
\item{peaks}{The user validated location (ratio Test/Normal) of the peaks before normalization.}
\item{ploidy}{The user validated ploidy of the peaks before normalization.}
        
}

\value{An object of class \code{"CNAnorm"}}

\author{Stefano Berri \email{s.berri@leeds.ac.uk}}
\examples{
data(CN)
# see peakPloidy documentation to know how object CN is created
CN <- validation(CN)
}
\seealso{
\code{\link{ratio.s}}, \code{\link{CNAnorm-class}}
}

\note{It is implicitly called by \code{discreteNorm} if no validation was manually performed}
