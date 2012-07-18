\name{gPar}
\docType{data}
\alias{gPar}
\title{An object with the default graphical parameters}

\description{
This data object is used by some plotting methods and contains the default 
values. User can change graphical parameters by changing this object

The object consists of several layers refering to different plots and different
properties. Here an indicative description

gPar$genome: prameters here refer to the plot produced by plotGenome

graphical parameters: see ?par

$colors: specify colors
$cex: specify relative size  - points, text...
$lwd: specify line width
$lty: specify line type - solid, dashed
$mar: specify margins



}
\usage{data(gPar)}
\format{A S3 object}
% \references{
% McNeil, D. R. (1977) \emph{Interactive Data Analysis}.
% New York: Wiley.
% }
\keyword{datasets, object}
