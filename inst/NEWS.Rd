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
    \item \code{profileMark()} was not outputting \code{crti.val} column
      when \code{crit.val = "bonf"} argument was chosen. This has now
      been amended.
      \item \code{quickEst()} has been improved and now provides an even
      more accurate and efficient distance calculations for linkage maps.
    }
}
}

