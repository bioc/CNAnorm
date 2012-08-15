\name{chromosomesPosition}
\docType{methods}
\alias{chrs}
\alias{chrs-methods}
\alias{chrs<-}
\alias{chrs<--methods}
\alias{arms}
\alias{arms-methods}
\alias{arms,CNAnorm-method}
\alias{pos}
\alias{pos-methods}
\alias{pos,CNAnorm-method}
\alias{pos<-}
\alias{pos<--methods}
\title{ Accessors methods for Function ratio in Package `CNAnorm' }
\description{
\code{chrs} returns/set the name of chromosomes/contigs

\code{arms} retruns the name of the chromosome and arm. A data frame containing
ideogram information has to be provided. See ?hg19_hs_ideogr for an example

\code{pos} returns/set the position of starting window. \bold{Be careful!} 
If you need to change data, it is better to change the input data and start over.
}

\usage{

chrs(object)
\S4method{pos}{CNAnorm}(object, show = "start")
\S4method{arms}{CNAnorm}(object, banding_df)
}

\arguments{
 
\item{object}{An object of Class \code{"CNAnorm"}}
\item{show}{The position to show: 'start', 'end'}
\item{banding_df}{A data frame with infromation about ideogram }
}

\value{\code{chrs} and \code{arms} return a character vector, \code{pos} returns a numeric vector}
\author{Stefano Berri <s.berri@leeds.ac.uk>}
\examples{
data(LS041)
data(hg19_hs_ideogr)
CN <- dataFrame2object(LS041)
dataFrameNames <- as.character(LS041$Chr)
objectNames <- chrs(CN)
armNames <- arms(CN, hg19_hs_ideogr)
# check the names are, indeed, the same
all(dataFrameNames == objectNames)
# make shorter names, drop the first three letters ('chr')
shortNames <- substr(chrs(CN),4,nchar(chrs(CN)))
chrs(CN) <- shortNames


# retrieve all new names
unique(chrs(CN))
unique(armNames)
}
\seealso{
\code{\link{gcNorm}}, \code{\link{CNAnorm-class}}, \code{\link{hg19_hs_ideogr}}
}



\keyword{methods, accessor, get}
