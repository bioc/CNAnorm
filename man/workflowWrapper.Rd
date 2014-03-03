\name{CNAnormWorkflow}
\alias{CNAnormWorkflow}
\title{Wrapper to \code{"CNAnorm"} workflow}
\description{This function is a wrapper to use for a fully automated CNAnorm workflow where
interactivity is not required. It contains MOST possible paramenters. Defaults
are set to to run a standard and conservative workflow.
}

\usage{
CNAnormWorkflow(dataFrame, gc.do=FALSE, gc.exclude=character(0),
    gc.maxNumPoints=10000, smooth.do=TRUE, smooth.lambda=7,
    smooth.other=list(), peak.method='closest', peak.exclude=character(0),
    peak.ploidyToTest=12, peak.sd=5, peak.dThresh=0.01, peak.n=2048,
    peak.adjust=.9, peak.force.smooth=TRUE, peak.reg=FALSE, peak.ds=1.5,
    peak.zero.count=FALSE, peak.other=list(), DNAcopy.do=TRUE,
    DNAcopy.independent.arms=FALSE, DNAcopy.ideograms=NULL,
    DNAcopy.smooth=list(), DNAcopy.segment=list(), DNAcopy.weight=character(),
    dNorm.normBy=NULL)

}

\arguments{
All arguments are explained in the relative functions
\item{dataFrame}{A data frame with columns Chr, Pos, Test, Norm and
    optional GC. See dataFrame2object}
\item{gc.do}{Specify if GC correction need to be done. See gcNorm }    
\item{gc.exclude}{See gcNorm}    
\item{gc.maxNumPoints}{See gcNorm}    
\item{smooth.do}{Specify if smoothing need to be done. See addSmooth}    
\item{smooth.lambda}{See addSmooth}    
\item{smooth.other}{A list of other parameters to pass to the smoothing
    function. See addSmooth}    
\item{peak.method}{See peakPloidy}    
\item{peak.exclude}{See peakPloidy}    
\item{peak.ploidyToTest}{See peakPloidy}    
\item{peak.sd}{See peakPloidy}    
\item{peak.dThresh}{See peakPloidy}    
\item{peak.n}{See peakPloidy}    
\item{peak.adjust}{See peakPloidy}    
\item{peak.force.smooth}{See peakPloidy}    
\item{peak.reg}{See peakPloidy}    
\item{peak.ds}{See peakPloidy}    
\item{peak.zero.count}{See peakPloidy}    
\item{peak.other}{A list of other parameters to be passed to funtions for peak
    detection See peakPloidy}    
\item{DNAcopy.do}{Specify if segmentation with DNAcopy need to be done. See addDNACopy}    
\item{DNAcopy.independent.arms}{See addDNACopy}    
\item{DNAcopy.ideograms}{See addDNACopy}    
\item{DNAcopy.smooth}{See addDNACopy}    
\item{DNAcopy.segment}{See addDNACopy}    
\item{DNAcopy.weight}{See addDNACopy}
\item{dNorm.normBy}{See discreteNorm}    
    
}

\value{
An object of Class \code{"CNAnorm"}
}

\author{Stefano Berri <s.berri@leeds.ac.uk>}

\examples{
data(LS041)
CN <- CNAnormWorkflow(LS041)
}

\seealso{
\code{\link{dataFrame2object}}, \code{\link{gcNorm}}, \code{\link{addSmooth}}, 
    \code{\link{peakPloidy}},  \code{\link{addDNACopy}}, \code{\link{discreteNorm}}
}

\keyword{workflow}
