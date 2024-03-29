% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{hitlist_elutriation}
\alias{hitlist_elutriation}
\title{Hitlist from the elutriation dataset}
\format{
## `hitlist_elutriation`
A data frame with 813 rows and 3 columns:
\describe{
  \item{id}{Uniprot protein ID}
  \item{treatment}{The treatment in which it was identified as a hit}
  \item{category}{Category of the hit}
}
}
\source{
This data were obtained from 7.	Dai, L., T. Zhao, X. Bisteau et al. (2018)
        Modulation of Protein-Interaction States through the Cell Cycle.
        Cell, 173, 1481-1494.e13.
}
\usage{
hitlist_elutriation
}
\description{
An hitlist from the elutriation dataset obtained with the \code{imprints_score} function
from the package \code{IMPRINTS.CETSA} using default parameters except useMAD set to TRUE.
}
\keyword{datasets}
