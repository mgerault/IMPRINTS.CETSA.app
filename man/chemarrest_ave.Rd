% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{chemarrest_ave}
\alias{chemarrest_ave}
\title{Chemarrest average fold-change dataset}
\format{
## `chemarrest_ave`
A data frame with 2777 rows and 23 columns:
\describe{
  \item{id}{Uniprot protein ID}
  \item{description}{Protein description}
  \item{37C_G1S}{Log2 fold-change G1S vs G1S phase at 37°C}
  \item{37C_PM}{Log2 fold-change PM vs G1S phase at 37°C}
  \item{37C_S}{Log2 fold-change S vs G1S phase at 37°C}
  \item{47C_G1S}{Log2 fold-change G1S vs G1S phase at 47°C}
  \item{47C_PM}{Log2 fold-change PM vs G1S phase at 47°C}
  \item{47C_S}{Log2 fold-change S vs G1S phase at 47°C}
  \item{50C_G1S}{Log2 fold-change G1S vs G1S phase at 50°C}
  \item{50C_PM}{Log2 fold-change PM vs G1S phase at 50°C}
  \item{50C_S}{Log2 fold-change S vs G1S phase at 50°C}
  \item{52C_G1S}{Log2 fold-change G1S vs G1S phase at 52°C}
  \item{52C_PM}{Log2 fold-change PM vs G1S phase at 52°C}
  \item{52C_S}{Log2 fold-change S vs G1S phase at 52°C}
  \item{54C_G1S}{Log2 fold-change G1S vs G1S phase at 54°C}
  \item{54C_PM}{Log2 fold-change PM vs G1S phase at 54°C}
  \item{54C_S}{Log2 fold-change S vs G1S phase at 54°C}
  \item{57C_G1S}{Log2 fold-change G1S vs G1S phase at 57°C}
  \item{57C_PM}{Log2 fold-change PM vs G1S phase at 57°C}
  \item{57C_S}{Log2 fold-change S vs G1S phase at 57°C}
  \item{sumUniPeps}{Median of the number of unique peptides from all temperatures}
  \item{sumPSMs}{Median of the PSMs from all temperatures}
  \item{countNum}{Median of the abundance count from all reference channel}
}
}
\source{
This data were obtained from 7.	Dai, L., T. Zhao, X. Bisteau et al. (2018)
        Modulation of Protein-Interaction States through the Cell Cycle.
        Cell, 173, 1481-1494.e13.
}
\usage{
chemarrest_ave
}
\description{
An IMPRINTS-CETSA dataset containing the mean log2 fold-change between different cell cycle phase.
}
\keyword{datasets}
