\name{exportTable}
\docType{methods}
\alias{exportTable}
\alias{exportTable-methods}
\alias{exportTable,CNAnorm-method}
\title{ Methods for Function exportTable in Package `CNAnorm' }
\description{
\code{exportTable} write a table with normalised values of each window. A wrapper to \code{"write.table"} 
}


\usage{
\S4method{exportTable}{CNAnorm}(object, file = "CNAnorm_table.tab", show = 'ratio', 
    sep = "\t", row.names = FALSE, ...)
}

\arguments{
 
\item{object}{an object of Class \code{"CNAnorm"}}
\item{file}{name of the file to save to}
\item{show}{what should be reported in the table: \code{"ratio"}: the normalized 
ratio (a value of 1 means diploid). \code{"ploidy"}: the same as ratio * 2.
\code{"center"}: report ratio centered on the most abbundant copy. Ratio of 1
means that the most abbundant ``state'' is centered to 1}
\item{sep}{the field separator string.}
\item{row.names}{either a logical value indicating whether the row number 
    should be written or a character vector of row names to be written.}
\item{...}{Extra arguments to be passed to \code{"write.table"}}
}

\value{An object of class \code{"CNAnorm"}}

\details{It produces a tab delimited text file with the following columns:

    Chr: Chromosome/contig name.

    Pos: Starting position of the window.
    
    Ratio: Ratio Test/Normal for each window after GC correction.
    
    Ratio.n: Ratio Test/Normal or ploidy for each window after normalisation.
    
    Ratio.s.n: Smoothed and normalised ratio Test/Normal or ploidy for each window.
    
    SegMean: Mean of the segment this window belongs to.
    
    SegMean.n: Normalised mean ratio Test/Normal or ploidy of the segment this window belongs to.


}
    

\author{Stefano Berri \email{s.berri@leeds.ac.uk}} 

\examples{
data(CN)
CN <- validation(CN)
CN <- discreteNorm(CN)
exportTable(CN, file = "CNAnorm_table.tab", show = 'ploidy')
}
\seealso{
    \code{\link{write.table}}
}



\keyword{methods, export, dataframe, table}
