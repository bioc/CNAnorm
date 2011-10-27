\name{retrieve peaks and ploidy}
\docType{methods}
\alias{sugg.peaks}
\alias{sugg.peaks-methods}
\alias{sugg.ploidy}
\alias{sugg.ploidy-methods}
\alias{valid.peaks}
\alias{valid.peaks-methods}
\alias{valid.ploidy}
\alias{valid.ploidy-methods}
\title{ Methods for Function to retrieve suggested/validated ploidy and peaks in Package `CNAnorm' }
\description{
\code{sugg.peaks} returns the location of peaks before normalization as 
    suggested by \code{peakPloidy}.

\code{sugg.ploidy} returns the ploidy of the peaks as suggested by
    \code{peakPloidy}.

\code{valid.peaks} returns the location of peaks before normalization as 
    validated after calling method \code{"validation"}
   
\code{valid.ploidy} returns the validated ploidy of the peaks as validated
    after calling method \code{"validation"} 
}

\usage{
sugg.peaks(object)
sugg.ploidy(object)
valid.peaks(object)
valid.ploidy(object)
}

\arguments{
 
\item{object}{An object of Class \code{"CNAnorm"} after method
    \code{"peakPloidy"} was called} }

\value{A numeric vector}
\author{Stefano Berri <s.berri@leeds.ac.uk>}
\examples{
data(CN)
# see peakPloidy documentation to know how object CN is created
plot(sugg.ploidy(CN), sugg.peaks(CN))

}
\seealso{
\code{\link{gcNorm}}, \code{\link{CNAnorm-class}}, \code{DNAcopy}
}



\keyword{methods, accessor, get, set}
