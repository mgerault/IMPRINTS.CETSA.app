% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/imprints_barplotting_peptides.R
\name{imprints_barplotting_peptides}
\alias{imprints_barplotting_peptides}
\title{imprints_barplotting_peptides}
\usage{
imprints_barplotting_peptides(
  data,
  treatmentlevel = get_treat_level(data),
  format = c("individual_peptide", "RESP_peptide", "peptide_one"),
  printBothName = TRUE,
  printGeneName = FALSE,
  witherrorbar = TRUE,
  withpoint = FALSE,
  pointperrep = TRUE,
  layout = NULL,
  colorpanel = PaletteWithoutGrey(treatmentlevel),
  usegradient = FALSE,
  colorgradient = c("#4575B4", "ivory", "#D73027"),
  linegraph = FALSE,
  log2scale = TRUE,
  ratio = 1.2,
  ret_plot = TRUE,
  save_pdf = FALSE,
  toplabel = "IMPRINTS-CETSA bar plotting",
  leftlabel = "",
  bottomlabel = "",
  pdfname = "barplot",
  pdfheight = 12,
  pdfwidth = 12
)
}
\arguments{
\item{data}{dataset after \code{imprints_sequence_peptides} to plot.}

\item{treatmentlevel}{a vector of treatment labels, such as c("DMSO","TNFa","AT26533")
the order determines the arrangement, so in this case DMSO
group would be the first group}

\item{format}{Format of the plot; either \code{individual_peptide} which will plot each peptide individually as
it would be done for each protein with \code{\link{imprints_barplotting_app}}, \code{RESP_peptide} which will plot
the RESP plot (you'll need exactly two 'peptides' per protein) or \code{peptide_one} which will plot all peptides
in one plot for each protein.
Default is \code{individual_peptide}.}

\item{printBothName}{A logical to tell if you want to print the both protein names on the plot}

\item{printGeneName}{A logical to tell if you want to print the gene names on the plot}

\item{witherrorbar}{A logical to print or not the error bar on the plot}

\item{withpoint}{A logical to print or not the data point of each replicate on the plot on top of the bars}

\item{pointperrep}{A logical to separate the point per replicate; only active when withpoint is set to TRUE}

\item{layout}{a vector indicating the panel layout for multi-panel plots per page,
default value is c(2,3) for set containing data, otherwise c(4,3), use when save_pdf = TRUE}

\item{colorpanel}{a vector of color scheme provided by default with the function PaletteWithoutGrey}

\item{usegradient}{whether the barplot should be draw in color gradient format}

\item{colorgradient}{the color scheme of gradient applied, default value c("#4575B4","ivory", "#D73027")}

\item{linegraph}{whether to plot the graph in a line graph format, default set to FALSE}

\item{log2scale}{whether the yscales should be in log2 scale, default set to TRUE}

\item{ratio}{aspect ratio of the plot, default set to 0.6}

\item{ret_plot}{Logical to tell if you want to return the last plot}

\item{save_pdf}{A logical to tell if you want to save plots in a pdf file}

\item{toplabel}{textual label at the top part of the page}

\item{leftlabel}{textual label at the left side of the page}

\item{bottomlabel}{textual label at the bottom part of the page}

\item{pdfname}{textual label of the pdf file}

\item{pdfheight}{a number indicate the height of pdf file, default value 12}

\item{pdfwidth}{a number indicate the width of pdf file, default value 12}
}
\value{
The imprints barplot
}
\description{
Function to generate IMPRINTS bar plot and pdf file with multipanel bar plots for IMPRINTS-CETSA data.
This function is based on the function imprints_barplotting from the IMPRINTS.CETSA package.
}
\seealso{
\code{\link{imprints_barplotting_app}}
}
