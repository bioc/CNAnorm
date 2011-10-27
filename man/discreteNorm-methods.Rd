\name{discreteNorm}
\docType{methods}
\alias{discreteNorm}
\alias{discreteNorm-methods}
\alias{discreteNorm,CNAnorm-method}
\title{ Methods for Function addSmooth in Package `CNAnorm' }
\description{
\code{discreteNorm} performs normalization of data using information on ploidy. 
    Implicitly calls \code{"validation"} if no validation was performed
}

\usage{
\S4method{discreteNorm}{CNAnorm}(object, normBy = object)
}

\arguments{
 
\item{object}{An object of Class \code{"CNAnorm"} to normalize}
\item{normBy}{An object of Class \code{"CNAnorm"} used to set normalization.
    It is possible, for instance, to normalize a sample at high resolution, using information
    obtained from the same sample at low resolution}
        
}

\references{ Gusnanto, A., Wood, H.M., Pawitan, Y., Rabbitts, P. and Berri, S.
(2011) \emph{Correcting for cancer genome size and tumour cell content enables
better estimation of copy number alterations from next generation sequence
data.} Bioinformatics}


\value{An object of class \code{"CNAnorm"}}

\author{Stefano Berri \email{s.berri@leeds.ac.uk} and 
    Arief Gusnanto \email{a.gusnanto@leeds.ac.uk}}
\examples{
data(CN)
# see peakPloidy documentation to know how object CN is created
CN <- discreteNorm(CN)

}
\seealso{
\code{\link{validation}}, \code{\link{peakPloidy}}
}

