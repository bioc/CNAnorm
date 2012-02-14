\name{peakPloidy}
\docType{methods}
\alias{peakPloidy}
\alias{peakPloidy-methods}
\alias{peakPloidy,CNAnorm-method}
\title{ Methods for Function peakPloidy in Package `CNAnorm' }
\description{
\code{peakPloidy} Estimate most likely ploidy of genome looking at distribution of smoothed ratio. 
}

\usage{
\S4method{peakPloidy}{CNAnorm}(object, method = 'mixture', exclude = character(0), 
    ploidyToTest = 12, sd = 5, dThresh = 0.01, n = 2048, adjust = .9, force.smooth = TRUE,
    reg = FALSE, ds = 1.5, zero.cont = FALSE, ...)
}

\arguments{
 
\item{object}{An object of Class \code{"CNAnorm"}}
\item{exclude}{A character vector with names of Chromosomes/Contigs
    \bold{not} to use to estimate ploidy.}
\item{method}{A character element matching either \code{"mixture"}, \code{"density"},
    \code{"median"}, \code{"mode"} or \code{"closest"}. \code{"mixture"} will
    fit a mixture model to find peaks; \code{"density"} will use the density
    function to find peaks; \code{"median"} \code{"mode"} and \code{"closest"}
    will only find one peak at the median, mode or peak closest to the median
    respectively. No tumour content or reliable estimated ploidy will be
    provided. These methods are ment to perform a more ``standard''
    normalisation, without stratching the data. Suggested for germline CNV 
    or a fully automated process that does not requires a normalisation on
    integer copy number or for highly heterogeneous sample where such normalisation
    would not be possible.
    Non ambigous partial matches can be used.}
\item{ploidyToTest}{Maximum ploidy allowed. \bold{Warnings!} Computation time
    increases exponentially with this parameter if using \code{"density"}.}
\item{adjust}{The \code{"adjust"} parameter passed to the \code{density} function.}
\item{n}{The \code{"n"} parameter passed to the \code{density} function.}
\item{force.smooth}{If the input object does not have smoothed ratio, it will smooth
     using \code{"addSmooth"}. It is highly recomended to use \code{"force.smooth = TRUE"}}
\item{sd}{Parameter to filter outliers. Values greater than i
    median + sd * standard deivations will be ignored while detecting peaks and ploidy.}
\item{dThresh}{Parameter to filter outliers. Values with a density lower than
    max(density)*dThresh will be ignored while detecting peaks and ploidy.}
\item{reg}{Parameter for mixture model:  If set TRUE, the starting point for EM
    algorithm will be optimized through several methods including regular grid on
    the ratio distribution. The default is FALSE, by which the starting values are
    taken from the quantiles of the distribution.}
\item{ds}{Parameter for mixture model: A constraint in EM algorithm of minimum
    distance between mean estimates, in terms of median standard deviation of the
    mixture components.} 
\item{zero.cont}{Parameter for mixture model: An argument for the mixture
    model. If set TRUE, the EM algorithm considers exactly-zero ratios as a mixture
    component.} 
\item{...}{Extra parameters to be passed to funtions for peak detection,
    specific to each of the methods (deinsity or mixture), se below for
    details.}

}

\note{
Other optional parameters to be passed (...)

\bold{mixture method}

\bold{density method}

\bold{peakRatio}{Threshold used to call a peak. Peaks smaller than 
    \code{maxPeakHight/peakRatio} are not considered peaks.}

\bold{spacingTolerance}{A parameter smaller than 1 which describes how
strict the program should be on alternative solutions.
Only solution with the best R^2 >= max(R^2)*spacingTolerance will be
considered as valid.}

\bold{interceptRatio}{Minimum value for the intercept of the linear model.
Ideally, should be zero, but the default allows a little flexibility.}

}

\value{An object of class \code{"CNAnorm"}}

\references{ Gusnanto, A., Wood, H.M., Pawitan, Y., Rabbitts, P. and Berri, S.
(2011) \emph{Correcting for cancer genome size and tumour cell content enables
better estimation of copy number alterations from next generation sequence
data.} Bioinformatics}

\author{Stefano Berri \email{s.berri@leeds.ac.uk} and Arief Gusnanto 
    \email{a.gusnanto@leeds.ac.uk}
}


\examples{
data(LS041)
CN <- dataFrame2object(LS041)
chr2skip <- c("chrY", "chrM")
CN <- gcNorm(CN, exclude = chr2skip)
CN <- addSmooth(CN, lambda = 3)
CN <- peakPloidy(CN, exclude = chr2skip)
# this object CN is what you obtain when you load 
# data(CN)
}


\seealso{
\code{\link{CNAnorm-class}}, \code{\link{density}}
}

\keyword{methods, normalization}
