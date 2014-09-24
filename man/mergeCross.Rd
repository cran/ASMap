\name{mergeCross}
\alias{mergeCross}
\title{
Merge linkage groups of an R/qtl cross object
}
\description{
Merges linkage groups of an R/qtl cross object from a user specified list.
}
\usage{
mergeCross(cross, merge = NULL, gap = 5)
}
\arguments{
  \item{cross}{
An R/qtl \code{cross} object with any class structure.
}
\item{merge}{
  A list with elements containing the linkage groups to be merged with
  each element named by the proposed linkage group name (see Examples).
}
\item{gap}{
The cM gap to put between the merged map elements in the complete linkage group.
}
}
\details{
This merging function allows you to perform multiple merges of two or more linkage groups
in one call. User should ensure linkage group names are correct and
that proposed linkage group names do not already exist.
}
\value{
The cross object is returned with identical class structure as the
inputted cross object. The \code{"geno"} element should now contain
merged linkage groups for the user defined merges.
}
\author{
Julian Taylor
}
\seealso{
\code{\link{breakCross}}
}
\examples{

data(mapDH, package = "ASMap")

mapDH1 <- breakCross(mapDH, split = list("4A" = "4A.m.8"))
pull.map(mapDH1)[["4A.1"]]
pull.map(mapDH1)[["4A.2"]]

mapDH2 <- mergeCross(mapDH1, merge = list("4A" = c("4A.1","4A.2")))
pull.map(mapDH2)[["4A"]]
}
\keyword{misc}