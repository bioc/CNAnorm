\name{addDNACopy}
\docType{methods}
\alias{addDNACopy}
\alias{addDNACopy-methods}
\alias{addDNACopy,CNAnorm-method}
\title{ Methods for Function addDNACopy in Package `CNAnorm' }
\description{
\code{addSmooth} segment ratio values in Package `CNAnorm' using DNACopy 
}
\section{Methods}{
\describe{

\item{\code{signature(object = "CNAnorm")}}{ Segment ratio values
    on an object of class \code{"CNAnorm"}. Returns the same object with 
    extra slots (\code{segMean}, \code{segID} )}
}
}

\usage{
\S4method{addDNACopy}{CNAnorm}(object, independent.arms = FALSE, ideograms = NULL)
}

\arguments{
 
\item{object}{An object of Class \code{"CNAnorm"}}
\item{independent.arms}{Boolean. If TRUE chromosomes arms will be treated as independent, ideograms must be provided}
\item{ideograms}{A data frame containing information about banding. See ?hg19_hs_ideogr for more information}
        
}

\value{An object of class \code{"CNAnorm"}}

\author{Stefano Berri \email{s.berri@leeds.ac.uk} and 
    Arief Gusnanto \email{a.gusnanto@leeds.ac.uk}}
\examples{
data(LS041)
CN <- dataFrame2object(LS041)
CN <- addDNACopy(CN)
}
\seealso{
\code{\link{segMean}}, \code{\link{CNAnorm-class}}, \code{DNAcopy}, \code{\link{hg19_hs_ideogr}}
}

\references{
Venkatraman, E. S. and Olshen, A. B. (2007) \emph{A faster circular binary
segmentation algorithm for the analysis of array CGH data}. Bioinformatics 
}


\keyword{methods, segmentation, smoothing}
