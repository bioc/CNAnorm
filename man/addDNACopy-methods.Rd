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
\S4method{addDNACopy}{CNAnorm}(object, independent.arms = FALSE, ideograms = NULL, 
    DNAcopy.smooth = list(), DNAcopy.weight = character(), DNAcopy.segment = list())
}

\arguments{
 
\item{object}{An object of Class \code{"CNAnorm"}}
\item{independent.arms}{Boolean. If TRUE chromosomes arms will be treated as independent, ideograms must be provided}
\item{ideograms}{A data frame containing information about banding. See ?hg19_hs_ideogr for more information}
\item{DNAcopy.smooth}{A list of parameters to be passed to function `smooth.CNA` in package DNAcopy}
\item{DNAcopy.weight}{A character value of one of these values. 'poisson',
    'gaussian', 'twoTailQuantile', 'oneTailQuantile'. It specifies a way to give
    weight to each window depending on how much coverage in the normal deviates
    from the median for that chromosome. Options are listed in decresing order
    of stringensy. See Details}
\item{DNAcopy.segment}{A list of parameters to be passed to function `segment` in package DNAcopy. Parameters `weights` and `verbose` are not accepted}
        
}

\details{
    'poisson': windows with coverage more or less than 2*sqrt(mean) from the
    mean are weighted down. The most stringent. Recommended for unbiased genome
    wide sequencing.

    'gaussian': windows with coverage more or less than 2*sd from the median
    are weighted down. Recommended for genome wide sequencing where coverage
    in the normal is far from poisson distribution.

    'twoTailQuantile': windows with coverage outside 5-95th quantile are weighted
    down. Recommended when coverage is far from a normal distribution - such as
    capture experiments -

    'oneTailQuantile': windows with coverage lower than 5th quantile are weighted
    down. Recommended when coverage is far from a normal distribution - such as
    capture experiments. Does not weight down windows with very high coverage.
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
