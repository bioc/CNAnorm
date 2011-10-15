\name{ratio}
\docType{methods}
\alias{ratio}
\alias{ratio-methods}
\alias{ratio.n}
\alias{ratio.n-methods}
\alias{ratio.s}
\alias{ratio.s-methods}
\alias{ratio.s.n}
\alias{ratio.s.n-methods}
\alias{segMean}
\alias{segMean-methods}
\alias{segMean.n}
\alias{segMean.n-methods}
\title{ Methods for Function ratio in Package `CNAnorm' }
\description{
\code{ratio} returns the Test/Normal ratio from an object of class \code{CNAnorm}.
    ratio is corrected for GC content if \code{gcNorm} was called.

\code{ratio.n} returns the Test/Normal \bold{normalized} ratio from an object of class \code{CNAnorm}
    after normalization. Its input is \code{ratio(object)}

\code{ratio.s} returns the Test/Normal \bold{smoothed} ratio from an object of class \code{CNAnorm}
    Its input is \code{ratio(object)}

\code{ratio.s.n} returns the Test/Normal \bold{smoothed and normalized} ratio 
    from an object of class \code{CNAnorm}. Its input is \code{ratio.s(object)}

\code{segMean} returns the mean of the segments as produced by DNACopy 

\code{segMean.n} returns the \bold{normalized} mean of the segments 

}

\usage{
ratio(object)
ratio.n(object)
ratio.s(object)
ratio.s.n(object)
segMean(object)
segMean.n(object)
}

\arguments{
 
\item{object}{An object of Class \code{"CNAnorm"}}
}

\value{A numeric vector}
\author{Stefano Berri <s.berri@leeds.ac.uk>}
\examples{
data(LS041)
CN <- dataFrame2object(LS041)
ratio.original <- ratio(CN)
CN.gcNorm <- gcNorm(CN, exclude = c("chrX", "chrY", "chrM"))
ratio.gc.corrected <- ratio(CN.gcNorm)

}
\seealso{
\code{\link{gcNorm}}, \code{\link{CNAnorm-class}}, \code{DNAcopy}
}



\keyword{methods, accessor, get, set}
