% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/imprints_remove_peptides.R
\name{imprints_remove_peptides}
\alias{imprints_remove_peptides}
\title{imprints_remove_peptides}
\usage{
imprints_remove_peptides(
  data,
  proteins = NULL,
  sequence,
  margin = 2,
  mode = c("remove", "keep")
)
}
\arguments{
\item{data}{A peptides dataset, typically after \code{imprints_sequence_peptides}.
It needs to contains the columns 'Master Protein Accessions' and 'Positions in Master Proteins'.}

\item{proteins}{The proteins from which you want to remove/keep the peptides.
If NULL and you put one sequence, it will remove/keep all these peptides for all proteins from peptides data.
If NULL and you put several sequence, the number of sequence you put needs to match the number
of proteins you have in your peptides data.}

\item{sequence}{The peptide position you want to remove/keep. If you put one sequence, it will remove/keep it for all proteins
from your peptides data; otherwise it needs to match the number of proteins you put or have in your dataset.
The format needs to be a number followed by a dash and another number, like this : '208-221'.}

\item{margin}{Numeric to tell the margin number of amino acid added to the sequences selected. This avoid
peptide selection ambigutiy when a peptide has only a small number of amino acid more than the sequence selected.
Default is set to 2.}

\item{mode}{Character to specify if you want to remove or only the peptides selected; either 'remove' or 'keep'.
Default is 'remove'.}
}
\value{
Your filtered dataset
}
\description{
Function to filter some peptides according there sequence. For example after using \code{imprints_sequence_peptides},
if you want to remove the peptides which corresponds to the cleaved sites.
}
