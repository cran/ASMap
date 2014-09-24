\name{mapDH}
\alias{mapDH}
\docType{data}
\title{A constructed linkage map for a doubled haploid wheat population}
\description{A constructed linkage map for a doubled
  haploid wheat population in the form of a constructed R/qtl object.
}
\usage{data(mapDH)}
\format{This data relates to a fully constructed linkage map of 599 markers genotyped on 218
  individuals. The linkage map consists of 23 linkage groups spanning
  the whole genome. Map distances were estimated using the Hidden Markov
  algorithm in \code{read.cross} with the \code{"kosambi"}
  mapping function. The map was originally constructed with
  MultiPoint and curated with MapManager and R/qtl. The data is in R/qtl
  format with a class structure \code{c("bc","cross")}. See
    \code{read.cross} documentation for more details on the
  format of this object.
}
\examples{

data(mapDH, package = "ASMap")

}
\keyword{datasets}