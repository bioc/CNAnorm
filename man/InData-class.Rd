\name{InData-class}
\Rdversion{1.1}
\docType{class}
\alias{InData-class}
\alias{length,InData-method}

\title{Class "InData" ~~~ }
\description{ A Class containing input data for CNA}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("InData", Chr, Pos, Test, Norm, ratio, GC)}.
}
\section{Slots}{
  \describe{
    \item{\code{Chr}:}{Object of class \code{"character"}. Name of the Chromosomes/Contigs of each window.}
    \item{\code{Pos}:}{Object of class \code{"numeric"}. Starting position of the each window.}
    \item{\code{Test}:}{Object of class \code{"numeric"}. Number of reads from Test in each window.}
    \item{\code{Norm}:}{Object of class \code{"numeric"}. Number of reads from Normal (Control) in each window.}
    \item{\code{ratio}:}{Object of class \code{"numeric"}. Ratio Test/Control
        in each window. Automatically computed if Test and Norm are provided, or
        user generated if Test and Norm are not know.}
    \item{\code{GC}:}{Object of class \code{"numeric"}. GC content of each window.}
  }
}
\section{Methods}{
  \describe{
    \item{length}{\code{signature(x = "InData")}: Returns number of windows from input data.}
	 }
}
\references{ Gusnanto, A., Wood, H.M., Pawitan, Y., Rabbitts, P. and Berri, S.
(2011) \emph{Correcting for cancer genome size and tumor cell content enables
better estimation of copy number alterations from next generation sequence
data.} Bioinformatics } 

\author{ Stefano Berri }


\seealso{
	\code{\linkS4class{CNAnorm}}
}
\examples{
data(LS041)
inObject <- new("InData", Chr = as.character(LS041$Chr), Pos = LS041$Pos,
    Test = LS041$Test, Norm = LS041$Norm, GC = LS041$GC)
CNA <- new("CNAnorm", InData = inObject) 
}
\keyword{classes}
