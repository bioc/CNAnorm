\name{dataFrame2object}
\alias{dataFrame2object}
\title{Convert a data frame into an object of Class \code{"CNAnorm"}}
\description{
Convert a data frame with column: Chr, Pos, Test, Norm and
    optional GC into object of class \code{"CNAnorm"}
}

\usage{
dataFrame2object(dataFrame)
}

\arguments{
\item{dataFrame}{A data frame with columns Chr, Pos, Test, Norm and
    optional GC}
}

\value{
An object of Class \code{"CNAnorm"}
}

\author{Stefano Berri <s.berri@leeds.ac.uk>}

\examples{
data(LS041)
CN <- dataFrame2object(LS041)
}

\seealso{
\code{\link{CNAnorm-class}}, \code{\link{InData-class}}, \code{\link{data.frame}}
}

\keyword{convert, data.frame}
