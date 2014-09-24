\name{ASMap-package}
\alias{ASMap-package}
\alias{ASMap}
\docType{package}
\title{
Additional functions for linkage map construction and manipulation of R/qtl
objects.
}
\description{
Additional functions for linkage map construction and manipulation of R/qtl
objects. This includes extremely fast linkage map clustering and marker
ordering using MSTmap (see Wu et al., 2008).
}
\details{
\tabular{ll}{
Package: \tab ASMap\cr
Type: \tab Package\cr
Version: \tab 0.3\cr
Date: \tab 2014-08-20\cr
License: \tab GPL 2\cr
}
Welcome to the first release version of the ASMap package!

One of the fundamental reasons why this package exists was to utilize
and implement the source code for the the Minimum Spanning Tree
algorithm derived in Wu et al. (2008) (reference 
below) for linkage map construction. The algorithm is lightning quick at
linkage group clustering and optimal marker ordering and can handle large numbers of
markers.

The package contains two very efficient functions, \code{mstmap.data.frame}
and \code{mstmap.cross}, that provide users with a highly flexible set
linkage map construction methods using the MSTmap
algorithm. \code{mstmap.data.frame} constructs a linkage map from
a data frame of genetic marker data and will use the entire
contents of the object to form linkage groups and optimally order
markers within each linkage group. \code{mstmap.cross} is a 
linkage map construction function for \pkg{qtl} package objects and can
be used to construct linkage maps in a flexible number of ways.
See \code{?mstmap.cross} for complete details.

To complement the computationally efficient linkage map construction
functions, the package also contains functions \code{pullCross} and
\code{pushCross} that allow the pulling/pushing markers of different
types to and from the linkage map. This system gives users the ability
to initially pull markers aside that are not needed for immediate
construction and push them back later if required. There are also
functions for fast numerical and graphical diagnosis of unconstructed
and constructed linkage maps. Specifically, there is an improved
\code{heatMap} that graphically displays pairwise recombination
fractions and LOD scores with separate legends for
each. \code{profileGen} can be used to simultaneously profile multiple
statistics such as recombination counts and double recombination
counts for individual lines across the constructed linkage
map. \code{profileMark} allows simultaneous graphical visualization of
marker or interval statistics profiles across the genome or subsetted for a predefined
set of linkage groups. Both of these functions utilize the power of
the advanced graphics package \pkg{lattice} to provide seamless multiple
displays. 

Other miscellaneous utilities for \pkg{qtl} objects include
\itemize{
  \item \code{mergeCross}: Merging of linkage groups
  \item \code{breakCross}: Breaking of linkage groups
  \item \code{combineMap}: Combining linkage maps
  \item \code{quickEst}: Very quick estimation of genetic map distances
  \item \code{genClones}: Reporting genotype clones
  \item \code{fixClones}: Consensus genotypes for clonal groups 
}

A comprehensive vignette showcasing the package will be available
shortly!
}
\author{
Julian Taylor, Dave Butler, Timothy Close, Yonghui Wu, Stefano Lonardi
Maintainer: Julian Taylor <julian.taylor@adelaide.edu.au>
}
\references{
Y. Wu, P. Bhat, T.J. Close, S. Lonardi, Efficient and Accurate
Construction of Genetic Linkage Maps from Minimum Spanning Tree of a
Graph Plos Genetics, Volume 4, Issue 10, 2008.
}
\keyword{package}
\seealso{
\code{\link[qtl]{qtl-package}}
}

