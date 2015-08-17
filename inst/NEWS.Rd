\name{news}
\title{Release information for ASMap}
\section{Changes in version 0.4-x}{
\subsection{README}{
  \itemize{
  \item This is the first documented NEWS release of ASMap. For this
  reason it contains information about features of the package as well
  as recent bug fixes.
 }
}
  \subsection{NEW FEATURES}{
    \itemize{
      \item The heat map in the package now uses a diverging colour
      palette chosen from the default palettes available in
      \pkg{RColorBrewer}.
      \item The package now includes a PDF vignette. This can be found
      on any operating system using the command

      \code{> vignette("ASMap")}

      The vignette includes extensive descriptions of the package
      functions as well as fully worked example that showcases the
      package functions throuh pre-diagnosis, construction and post
      diagnosis of a linkage map. It also describes how the package can
      be used efficiently for techniques such as fine mapping or combining
      maps. 
      \item The package now includes three fully documented marker data
      sets for users to explore. See the vignette for more details.    
}
}
\subsection{BUG FIXES}{
  \itemize{
     \item Fixed bug that returned inflated distances for \code{"ARIL"}
     populations when using \code{mstmap.data.frame()}.
     \item If the constructed map contained linkage groups with one
      marker \code{profileMark()} graphically displaced the interval
      statistics. This has now been amended.
     \item \code{pushCross()} can now push large volumes of co-located
      markers back into the constructed map very efficiently.   
    \item \code{profileMark()} was not outputting \code{crti.val} column
      when \code{crit.val = "bonf"} argument was chosen. This has now
      been amended.
      \item \code{quickEst()} has been improved and now provides an even
      more accurate and efficient distance calculations for linkage maps.
    }
}
}


