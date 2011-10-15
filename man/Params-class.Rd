\name{Params-class}
\Rdversion{1.1}
\docType{class}
\alias{Params-class}

\title{Class "Params"}
\description{ A Class containing some Parameters used in the analysis. Not heavily used at the moment}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Params")}, it is usually
    iniziated and populated with methods and functions of class \code{CNAnorm}
}  
\section{Slots}{
  \describe{
    \item{\code{method}:}{variable of class \code{"character"}. Record if the \code{peakPloidy}
    function was called using \code{density} or \code{mixture}}
    \item{\code{density.n}:}{The varaible \code{"n"} used when calling \code{peakPloidy}. 
        This variable is saved so that can be used later for drawing plots}
    \item{\code{density.adjust}:}{The varaible \code{"adjust"} used when calling \code{peakPloidy}. 
        This variable is saved so that can be used later for drawing plots}
    \item{\code{gc.excludeFromGCNorm}:}{Vector of class \code{"character"}.
    Name of the Chromosomes/Contigues not used for GC content correction.}
    \item{\code{gc.maxNumPoints}:}{One element vector of class \code{"numeric"}. 
    Specify how many points to use for GC correction }
    \item{\code{gp.excludeFromDensity}:}{Vector of class \code{"character"}. 
    Name of the Chromosomes/Contigues not used for peak guessing }
  }
}
\section{Methods}{
  \describe{
    \item{length}{\code{signature(x = "Params")}}
	 }
}
\references{ To be added }
\author{ Stefano Berri}


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
