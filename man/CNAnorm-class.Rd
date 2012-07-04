\name{CNAnorm-class}
\Rdversion{1.1}
\docType{class}
\alias{CNAnorm}
\alias{CNAnorm-class}
\alias{chrs,CNAnorm-method}
\alias{chrs<-,CNAnorm-method}
\alias{length,CNAnorm-method}
\alias{[,CNAnorm-method}
\alias{pos<-,CNAnorm-method}
\alias{length,CNAnorm-method}
\alias{ratio,CNAnorm-method}
\alias{ratio.n,CNAnorm-method}
\alias{ratio.s,CNAnorm-method}
\alias{ratio.s.n,CNAnorm-method}
\alias{segMean,CNAnorm-method}
\alias{segMean.n,CNAnorm-method}
\alias{sugg.peaks,CNAnorm-method}
\alias{sugg.ploidy,CNAnorm-method}
\alias{valid.peaks,CNAnorm-method}
\alias{valid.ploidy,CNAnorm-method}

\title{Class "CNAnorm"}
\description{Class to Contain and Describe copy number aberration (CNA) data from low coverage (approx 0.01 - 0.5X) Next Generation Sequencing}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CNAnorm", InData)}.
}

\section{Slots}{
  \describe{
    \item{\code{InData}:}{Object of class \code{"InData"}. Contains input data
        provided by the user. All slots have same length. Each element describe one window.
        See Class \code{"InData"}}
    \item{\code{DerivData}:}{Object of class \code{"DerivData"}. Contains data derived
        from \code{"InData"}. It can be Retrieved by the user, but methods should be used to populate
        \code{"DerivData"}. All slots have same length as input data. See Class DerivData}
    \item{\code{Res}:}{Object of class \code{"Res"}. Contains slots with obtained results. See Class \code{"Res"} }
    \item{\code{Params}:}{Object of class \code{"Params"}. Contains crucial parameters passed to some 
        of the methods for reusing in later steps or for documentation. }
  }
}
\section{Methods}{
  Summary of methods for class \code{"CNAnorm"}. Type \code{"methods ? methodName"} for more details about methodName.
  \describe{
    \item{chrs}{\code{signature(object = "CNAnorm")}: Retrieve Chromosomes/contig name }
    \item{chrs<-}{\code{signature(object = "CNAnorm")}: Set Chromosomes/contig name }
    \item{guessPeaksAndPloidy}{\code{signature(object = "CNAnorm")}: Estimate ploidy 
        of the sample, tumor content and other results saved in Slot \code{"Res"} }
    \item{length}{\code{signature(x = "CNAnorm")}: Returns number of element/windows }
    \item{[}{\code{signature(x = "CNAnorm")}: Produce on object of class \code{"CNAnorm"} with a subser of windows }
    \item{plotGenome}{\code{signature(object = "CNAnorm")}: Plot annotated normalized genome copy number }
    \item{plotPeaks}{\code{signature(object = "CNAnorm")}: Plot peaks and estimated/validated ploidy }
    \item{pos}{\code{signature(object = "CNAnorm")}: Retrieve Chromosomes/contig position }
    \item{pos<-}{\code{signature(object = "CNAnorm")}: Set Chromosomes/contig position }
    \item{ratio}{\code{signature(object = "CNAnorm")}: Retrieve ratio (Test/Control). 
        If gcNorm was called, ratio is GC normalized }
    \item{ratio.n}{\code{signature(object = "CNAnorm")}: Retrieve normalized ratio (not smoothed) }
    \item{ratio.s}{\code{signature(object = "CNAnorm")}: Retrieve smoothed ratio }
    \item{ratio.n.s}{\code{signature(object = "CNAnorm")}: Retrieve normalized smoothed ratio }
    \item{segMean}{\code{signature(object = "CNAnorm")}: Retrieve segmented ratio (as provided by DNAcopy) }
    \item{segMean.n}{\code{signature(object = "CNAnorm")}: Retrieve normalized segmented ratio }
	 }
}
\references{ CNA-norm: Discrete Normalization of Copy Number Alteration data from clinical samples (in preparation) }
\author{ Stefano Berri \email{s.berri@leeds.ac.uk} and Arief Gusnanto \email{a.gusnanto@leeds.ac.uk} }


\seealso{
	\code{\linkS4class{InData}}, \code{\linkS4class{DerivData}} for documentation on the slots.
}
\examples{
data(LS041)
CNA <- new("CNAnorm", InData = new("InData", Chr = as.character(LS041$Chr), Pos = LS041$Pos,
    Test = LS041$Test, Norm = LS041$Norm, GC = LS041$GC))
}
\keyword{classes}
