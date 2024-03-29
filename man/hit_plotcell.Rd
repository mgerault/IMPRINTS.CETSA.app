% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hitplot_cell.R
\name{hit_plotcell}
\alias{hit_plotcell}
\title{hit_plotcell}
\usage{
hit_plotcell(
  data,
  tit = "PI3K data in the cell",
  cond = NULL,
  cat_col_list = list(CC = "red", CN = "lightblue", NC = "yellow", ND = "#747474", NN =
    "#CCCCCC")
)
}
\arguments{
\item{data}{The output from hit_for_cell function}

\item{tit}{The title of the plot}

\item{cond}{The treatments you want to keep. If NULL, will take it all}

\item{cat_col_list}{A list containing the color for each category (CC, CN, NC, ND, NN)}
}
\value{
An interactive plot
}
\description{
Function to return the 'cell plot' of the proteins according to their category.
}
\seealso{
\code{\link{hit_for_cell}}
}
