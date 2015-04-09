\name{mapBCu}
\alias{mapBCu}
\docType{data}
\title{An unconstructed marker set for a backcross barley population}
\description{An unconstructed marker set for a backcross
  barley population in the form of an R/qtl object.
}
\usage{data(mapBCu)}
\format{This data relates to an unconstructed version of
  \code{\link{mapBC}} and consists of 3023 markers genotyped on 326
  individuals with markers randomly assorted on one large linkage group.
  The data is in R/qtl format with a class structure
  \code{c("bc","cross")}. See \code{read.cross()} documentation for more
  details on the format of this object. This data set forms the basis of
  the worked example in Chapter 3 of the vignette (see vignette("ASMap")
  for complete details
}
\examples{

data(mapBCu, package = "ASMap")

}
\keyword{datasets}