\name{addSmooth}
\docType{methods}
\alias{addSmooth}
\alias{addSmooth-methods}
\alias{addSmooth,CNAnorm-method}
\title{ Methods for Function addSmooth in Package `CNAnorm' }
\description{
\code{addSmooth} segment and smooth perform ratio values in Package `CNAnorm' 
}
\section{Methods}{
\describe{

\item{\code{signature(object = "CNAnorm")}}{ Segment and smooth perform ratio values
    on an object of class \code{"CNAnorm"}. Returns the same object with 
    extra slot (\code{ratio.s})}
}
}

\usage{
\S4method{addSmooth}{CNAnorm}(object, lambda = 7, ...)
}

\arguments{
 
\item{object}{An object of Class \code{"CNAnorm"}}
\item{lambda}{Degree of smoothness. See reference for more details}
\item{...}{Further arguments to pass to the function \code{smoothseg}}
        
}

\value{An object of class \code{"CNAnorm"}}

\author{Stefano Berri \email{s.berri@leeds.ac.uk} and 
    Arief Gusnanto \email{a.gusnanto@leeds.ac.uk}}
\examples{
data(LS041)
CN <- dataFrame2object(LS041)
CN.gcNorm <- gcNorm(CN, exclude = c("chrX", "chrY", "chrM"))
CN.smooth <- addSmooth(CN)
}
\seealso{
\code{\link{ratio.s}}, \code{\link{CNAnorm-class}}
}

\references{
Huang, J., Gusnanto, A., O'Sullivan, K., Staaf, J., Borg, A. and Pawitan, Y. 
(2007) \emph{Robust smooth segmentation approach for array CGH data analysis}.
Bioinformatics

}


\keyword{methods, segmentation, smoothing}
