\name{DerivData-class}
\Rdversion{1.1}
\docType{class}
\alias{DerivData-class}
\alias{length,DerivData-method}
\title{Class "DerivData"}
\description{ A Class containing data derived from InData used for further computation and plottting}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("DerivData")}, however DerivData is tipically 
populated using methods. If a slot has not been populated yet, it has zero length, otherwise slots
have the same length as \code{InData}
}
\section{Slots}{
  \describe{
    \item{\code{ratio}:}{Numeric vector with ratio Test/Normal. Optionally GC corrected}
    \item{\code{ratio.s}:}{Numeric vector with smoothed ratio.}
    \item{\code{ratio.n}:}{Numeric vector with normalized ratio.}
    \item{\code{ratio.s.n}:}{Numeric vector with normalized and smoothed ratio}
    \item{\code{segID}:}{Numeric vecotr with ID of segmented data (as provided by DNACopy).
        Each \code{segment} has a different ID.}
    \item{\code{segMean}:}{Numeric vector with mean value of the segment (as provided by DNACopy)}
    \item{\code{segMean.n}:}{Numeric vector with normalized segMean.}
    \item{\code{ok4density}:}{Logical vector. Specify wich values have been used to calculate density}
  }
}
\section{Methods}{
  \describe{
    \item{length}{\code{signature(x = "DerivData")}: Returns number of windows}
	 }
}
\references{ To be added }
\author{ Stefano Berri and Arief Gusnanto }


\seealso{
	\code{\linkS4class{CNAnorm}}, \code{\link{InData-class}}
}
\examples{
data(LS041)
inObject <- new("InData", Chr = as.character(LS041$Chr), Pos = LS041$Pos,
    Test = LS041$Test, Norm = LS041$Norm, GC = LS041$GC)
CNA <- new("CNAnorm", InData = inObject) 
}
\keyword{classes}
