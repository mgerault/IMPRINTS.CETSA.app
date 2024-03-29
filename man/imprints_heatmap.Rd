% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/imprints_heatmap.R
\name{imprints_heatmap}
\alias{imprints_heatmap}
\title{imprints_heatmap}
\usage{
imprints_heatmap(
  data,
  hit_summary = NULL,
  NN_data = NULL,
  PRcomplex_data = NULL,
  treatment,
  max_na = 0,
  response = c("both", "S", "D"),
  select_cat = c("CC", "CN", "NC"),
  saveHeat = FALSE,
  file_type = c("png", "pdf"),
  file_name = "Heatmap",
  titleH = "Elutriation heatmap",
  gradient_color = c("#005EFF", "#FFFFFF", "#FF0000"),
  cat_color = list(CC = "#FB4F0B", CN = "#0FAEB9", NC = "#E7B700"),
  back_color = "#FFFFFF",
  border_color = NULL
)
}
\arguments{
\item{data}{Output data from imprints_average}

\item{hit_summary}{The summary file from the hitlist output}

\item{NN_data}{The NN file from the hitlist output}

\item{PRcomplex_data}{Output data from imprints_complex_mapping. If not NULL and hit_summary NULL, will
print different heatmaps according to the protein complex.}

\item{treatment}{A character telling the condition from which you want to see the heatmap}

\item{max_na}{An integer indicating the maximum number of missing values per row (per protein)}

\item{response}{A character to tell if you want to the destabilized proteins, stabilized or both.
Accepted value are 'D' for destabilize, 'S' for stabilized or 'both'.}

\item{select_cat}{A character vector indicating the categories from which you want to see the heatmap}

\item{saveHeat}{Logical to tell if you want to save the heatmap}

\item{file_type}{The format file, either 'png' or 'pdf'.}

\item{file_name}{The file name}

\item{titleH}{The title for your heatmap}

\item{gradient_color}{The color for the gradient of the heatmap. Can only be of length three.}

\item{cat_color}{A list which contains the colors for each categories you selected.}

\item{back_color}{The color from the background of the heatmap (can be NULL)}

\item{border_color}{The color from the border of the plot (can be NULL)}
}
\value{
A grob object, the heatmap.
}
\description{
Function to get the heatmap from your data.
}
\seealso{
\code{\link{imprints_average}} , \code{\link{imprints_complex_mapping}}
}
