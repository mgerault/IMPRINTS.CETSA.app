% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PTMs_mapping.R
\name{imprints_ptms_peptides}
\alias{imprints_ptms_peptides}
\title{imprints_ptms_peptides}
\usage{
imprints_ptms_peptides(
  data,
  data_cleaved,
  control,
  minref = 2,
  PTM_FC_cutoff = 0.25,
  PTM_direction = c("hyper", "hypo", "both"),
  save_xlsx = TRUE,
  xlsxname = "RESP_PTMs_mapping"
)
}
\arguments{
\item{data}{The normalized peptides data set, i.e. the outpout from \code{imprints_normalize_peptides}.}

\item{data_cleaved}{The cleavage hits data set, i.e. the outpout from \code{imprints_cleaved_peptides}.}

\item{control}{The control treatment from your dataset.}

\item{minref}{The minimum number of references reporting the PTM in PhosphoSitePlus; Default is 2}

\item{PTM_FC_cutoff}{The minimum value for the maximum log2 fold-change of a peptide to be considered
significantly modified}

\item{PTM_direction}{Character to set if you look for PTMs in higher or lower abundance or both you can
select respectively hyper, hypo or both. Default is hyper.
Since the function assume you didn't add the PTM as a dynamic modification during your protein identification search,
if a peptide is "hypermodified' then peptide with negative value will be selected and conversely if a peptide is
'hypomodified' peptide with positive value will be kept.}

\item{save_xlsx}{Logical to tell if you want to save the categorized hits in an xlsx file.
Default to TRUE.}

\item{xlsxname}{The name of your saved file.}
}
\value{
A dataframe containing the potential proteins having PTMs
 instead of the protein being cleaved
}
\description{
Function to check if the hits found by \code{\link{imprints_cleaved_peptides}} could be due to PTMs and
 not  to protein cleavage.
}
\details{
When a hit is returned by \code{\link{imprints_cleaved_peptides}}, it means that this protein has a peptide position
where the IMPRINTS profiles of the two obtained parts are significantly different. This difference can be caused
by protein modification and mainly proteolysis but other more common PTMs can cause such a difference like
phosphorylation or ubiquitination. The aim here is to refilter the hit list and give the possible proteins which
could be modified based on the output of \code{\link{imprints_cleaved_peptides}}.
}
\seealso{
\code{\link{imprints_cleaved_peptides}}
}
