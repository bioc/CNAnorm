\name{gcNorm}
\docType{methods}
\alias{gcNorm}
\alias{gcNorm-methods}
\alias{gcNorm,CNAnorm-method}
\title{ Methods for Function gcNorm in Package `CNAnorm' }
\description{
\code{gcNorm} perform GC content normalization on ratio Test/Normal in Package `CNAnorm' 
}
\section{Methods}{
\describe{

\item{\code{signature(object = "CNAnorm")}}{ Perform GC content correction on 
    an object of class \code{"CNAnorm"}. Returns the same object with corrected ratio }
}
}

\usage{
\S4method{gcNorm}{CNAnorm}(object, exclude = character(0), maxNumPoints = 10000)
}

\arguments{
 
\item{object}{An object of Class \code{"CNAnorm"}}
\item{exclude}{A character vector with name of chromosomes/contigues 
    not to use to calculate GC content correction. All genome, however, will be corrected}
\item{maxNumPoints}{Maximum number of data points to fit the loess correction.
        For computational pourposes, if the number of points in \code{ratio(object)} is
        greater than \code{maxNumPoints}, only \code{maxNumPoints} randomly selected 
        will be used} 
}
\value{An object of class \code{"CNAnorm"}}
\author{Stefano Berri <s.berri@leeds.ac.uk>}
\examples{
data(LS041)
CN <- dataFrame2object(LS041)
# correct for GC content, but ignoring data from sex chromosomes and
# mitocondria
CN.gcNorm <- gcNorm(CN, exclude = c("chrX", "chrY", "chrM"))

}
\seealso{
\code{\link{loess}}, \code{\link{CNAnorm-class}}, \code{\link{ratio}}
}



\keyword{methods, normalization}
