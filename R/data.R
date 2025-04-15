#' Elutriation fold-change dataset
#'
#' An IMPRINTS-CETSA dataset containing the log2 fold-change between each bio replicate from
#' different cell cycle phase.
#'
#' @format ## `elutriation`
#' A data frame with 2731 rows and 59 columns:
#' \describe{
#'   \item{id}{Uniprot protein ID}
#'   \item{description}{Protein description}
#'   \item{37C_B1_G1}{Log2 fold-change G1 vs G1 phase at 37°C in the first bioreplicate}
#'   \item{37C_B1_G2}{Log2 fold-change G2 vs G1 phase at 37°C in the first bioreplicate}
#'   \item{37C_B1_S}{Log2 fold-change S vs G1 phase at 37°C in the first bioreplicate}
#'   \item{37C_B2_G1}{Log2 fold-change G1 vs G1 phase at 37°C in the second bioreplicate}
#'   \item{37C_B2_G2}{Log2 fold-change G2 vs G1 phase at 37°C in the second bioreplicate}
#'   \item{37C_B2_S}{Log2 fold-change S vs G1 phase at 37°C in the second bioreplicate}
#'   \item{37C_B3_G1}{Log2 fold-change G1 vs G1 phase at 37°C in the third bioreplicate}
#'   \item{37C_B3_G2}{Log2 fold-change G2 vs G1 phase at 37°C in the third bioreplicate}
#'   \item{37C_B3_S}{Log2 fold-change S vs G1 phase at 37°C in the third bioreplicate}
#'   \item{47C_B1_G1}{Log2 fold-change G1 vs G1 phase at 47°C in the first bioreplicate}
#'   \item{47C_B1_G2}{Log2 fold-change G2 vs G1 phase at 47°C in the first bioreplicate}
#'   \item{47C_B1_S}{Log2 fold-change S vs G1 phase at 47°C in the first bioreplicate}
#'   \item{47C_B2_G1}{Log2 fold-change G1 vs G1 phase at 47°C in the second bioreplicate}
#'   \item{47C_B2_G2}{Log2 fold-change G2 vs G1 phase at 47°C in the second bioreplicate}
#'   \item{47C_B2_S}{Log2 fold-change S vs G1 phase at 47°C in the second bioreplicate}
#'   \item{47C_B3_G1}{Log2 fold-change G1 vs G1 phase at 47°C in the third bioreplicate}
#'   \item{47C_B3_G2}{Log2 fold-change G2 vs G1 phase at 47°C in the third bioreplicate}
#'   \item{47C_B3_S}{Log2 fold-change S vs G1 phase at 47°C in the third bioreplicate}
#'   \item{50C_B1_G1}{Log2 fold-change G1 vs G1 phase at 50°C in the first bioreplicate}
#'   \item{50C_B1_G2}{Log2 fold-change G2 vs G1 phase at 50°C in the first bioreplicate}
#'   \item{50C_B1_S}{Log2 fold-change S vs G1 phase at 50°C in the first bioreplicate}
#'   \item{50C_B2_G1}{Log2 fold-change G1 vs G1 phase at 50°C in the second bioreplicate}
#'   \item{50C_B2_G2}{Log2 fold-change G2 vs G1 phase at 50°C in the second bioreplicate}
#'   \item{50C_B2_S}{Log2 fold-change S vs G1 phase at 50°C in the second bioreplicate}
#'   \item{50C_B3_G1}{Log2 fold-change G1 vs G1 phase at 50°C in the third bioreplicate}
#'   \item{50C_B3_G2}{Log2 fold-change G2 vs G1 phase at 50°C in the third bioreplicate}
#'   \item{50C_B3_S}{Log2 fold-change S vs G1 phase at 50°C in the third bioreplicate}
#'   \item{52C_B1_G1}{Log2 fold-change G1 vs G1 phase at 52°C in the first bioreplicate}
#'   \item{52C_B1_G2}{Log2 fold-change G2 vs G1 phase at 52°C in the first bioreplicate}
#'   \item{52C_B1_S}{Log2 fold-change S vs G1 phase at 52°C in the first bioreplicate}
#'   \item{52C_B2_G1}{Log2 fold-change G1 vs G1 phase at 52°C in the second bioreplicate}
#'   \item{52C_B2_G2}{Log2 fold-change G2 vs G1 phase at 52°C in the second bioreplicate}
#'   \item{52C_B2_S}{Log2 fold-change S vs G1 phase at 52°C in the second bioreplicate}
#'   \item{52C_B3_G1}{Log2 fold-change G1 vs G1 phase at 52°C in the third bioreplicate}
#'   \item{52C_B3_G2}{Log2 fold-change G2 vs G1 phase at 52°C in the third bioreplicate}
#'   \item{52C_B3_S}{Log2 fold-change S vs G1 phase at 52°C in the third bioreplicate}
#'   \item{54C_B1_G1}{Log2 fold-change G1 vs G1 phase at 54°C in the first bioreplicate}
#'   \item{54C_B1_G2}{Log2 fold-change G2 vs G1 phase at 54°C in the first bioreplicate}
#'   \item{54C_B1_S}{Log2 fold-change S vs G1 phase at 54°C in the first bioreplicate}
#'   \item{54C_B2_G1}{Log2 fold-change G1 vs G1 phase at 54°C in the second bioreplicate}
#'   \item{54C_B2_G2}{Log2 fold-change G2 vs G1 phase at 54°C in the second bioreplicate}
#'   \item{54C_B2_S}{Log2 fold-change S vs G1 phase at 54°C in the second bioreplicate}
#'   \item{54C_B3_G1}{Log2 fold-change G1 vs G1 phase at 54°C in the third bioreplicate}
#'   \item{54C_B3_G2}{Log2 fold-change G2 vs G1 phase at 54°C in the third bioreplicate}
#'   \item{54C_B3_S}{Log2 fold-change S vs G1 phase at 54°C in the third bioreplicate}
#'   \item{57C_B1_G1}{Log2 fold-change G1 vs G1 phase at 57°C in the first bioreplicate}
#'   \item{57C_B1_G2}{Log2 fold-change G2 vs G1 phase at 57°C in the first bioreplicate}
#'   \item{57C_B1_S}{Log2 fold-change S vs G1 phase at 57°C in the first bioreplicate}
#'   \item{57C_B2_G1}{Log2 fold-change G1 vs G1 phase at 57°C in the second bioreplicate}
#'   \item{57C_B2_G2}{Log2 fold-change G2 vs G1 phase at 57°C in the second bioreplicate}
#'   \item{57C_B2_S}{Log2 fold-change S vs G1 phase at 57°C in the second bioreplicate}
#'   \item{57C_B3_G1}{Log2 fold-change G1 vs G1 phase at 57°C in the third bioreplicate}
#'   \item{57C_B3_G2}{Log2 fold-change G2 vs G1 phase at 57°C in the third bioreplicate}
#'   \item{57C_B3_S}{Log2 fold-change S vs G1 phase at 57°C in the third bioreplicate}
#'   \item{sumUniPeps}{Median of the number of unique peptides from all temperatures}
#'   \item{sumPSMs}{Median of the PSMs from all temperatures}
#'   \item{countNum}{Median of the abundance count from all reference channel}
#' }
#' @source This data were obtained from 7.	Dai, L., T. Zhao, X. Bisteau et al. (2018)
#'         Modulation of Protein-Interaction States through the Cell Cycle.
#'         Cell, 173, 1481-1494.e13.
"elutriation"


#' Elutriation average fold-change dataset
#'
#' An IMPRINTS-CETSA dataset containing the mean log2 fold-change between different cell cycle phase.
#'
#' @format ## `elutriation_ave`
#' A data frame with 2731 rows and 23 columns:
#' \describe{
#'   \item{id}{Uniprot protein ID}
#'   \item{description}{Protein description}
#'   \item{37C_G1}{Log2 fold-change G1 vs G1 phase at 37°C}
#'   \item{37C_G2}{Log2 fold-change G2 vs G1 phase at 37°C}
#'   \item{37C_S}{Log2 fold-change S vs G1 phase at 37°C}
#'   \item{47C_G1}{Log2 fold-change G1 vs G1 phase at 47°C}
#'   \item{47C_G2}{Log2 fold-change G2 vs G1 phase at 47°C}
#'   \item{47C_S}{Log2 fold-change S vs G1 phase at 47°C}
#'   \item{50C_G1}{Log2 fold-change G1 vs G1 phase at 50°C}
#'   \item{50C_G2}{Log2 fold-change G2 vs G1 phase at 50°C}
#'   \item{50C_S}{Log2 fold-change S vs G1 phase at 50°C}
#'   \item{52C_G1}{Log2 fold-change G1 vs G1 phase at 52°C}
#'   \item{52C_G2}{Log2 fold-change G2 vs G1 phase at 52°C}
#'   \item{52C_S}{Log2 fold-change S vs G1 phase at 52°C}
#'   \item{54C_G1}{Log2 fold-change G1 vs G1 phase at 54°C}
#'   \item{54C_G2}{Log2 fold-change G2 vs G1 phase at 54°C}
#'   \item{54C_S}{Log2 fold-change S vs G1 phase at 54°C}
#'   \item{57C_G1}{Log2 fold-change G1 vs G1 phase at 57°C}
#'   \item{57C_G2}{Log2 fold-change G2 vs G1 phase at 57°C}
#'   \item{57C_S}{Log2 fold-change S vs G1 phase at 57°C}
#'   \item{sumUniPeps}{Median of the number of unique peptides from all temperatures}
#'   \item{sumPSMs}{Median of the PSMs from all temperatures}
#'   \item{countNum}{Median of the abundance count from all reference channel}
#' }
#' @source This data were obtained from 7.	Dai, L., T. Zhao, X. Bisteau et al. (2018)
#'         Modulation of Protein-Interaction States through the Cell Cycle.
#'         Cell, 173, 1481-1494.e13.
"elutriation_ave"


#' Hitlist from the elutriation dataset
#'
#' An hitlist from the elutriation dataset obtained with the \code{imprints_score} function
#' from the package \code{IMPRINTS.CETSA} using default parameters except useMAD set to TRUE.
#'
#' @format ## `hitlist_elutriation`
#' A data frame with 813 rows and 3 columns:
#' \describe{
#'   \item{id}{Uniprot protein ID}
#'   \item{treatment}{The treatment in which it was identified as a hit}
#'   \item{category}{Category of the hit}
#' }
#' @source This data were obtained from 7.	Dai, L., T. Zhao, X. Bisteau et al. (2018)
#'         Modulation of Protein-Interaction States through the Cell Cycle.
#'         Cell, 173, 1481-1494.e13.
"hitlist_elutriation"


#' The non hits from the elutriation dataset
#'
#' The non hits derived from the hitlist of the elutriation dataset
#'
#' @format ## `NN_elutriation`
#' A data frame with 4567 rows and 4 columns:
#' \describe{
#'   \item{id}{Uniprot protein ID}
#'   \item{description}{Protein description}
#'   \item{treatment}{The treatment in which it was identified as a hit}
#'   \item{category}{Category of the non hit, i.e. NN}
#' }
#' @source This data were obtained from 7.	Dai, L., T. Zhao, X. Bisteau et al. (2018)
#'         Modulation of Protein-Interaction States through the Cell Cycle.
#'         Cell, 173, 1481-1494.e13.
"NN_elutriation"




#' Chemarrest fold-change dataset
#'
#' An IMPRINTS-CETSA dataset containing the log2 fold-change between each bio replicate from
#' different cell cycle phase.
#'
#' @format ## `chemarrest`
#' A data frame with 2777 rows and 59 columns:
#' \describe{
#'   \item{id}{Uniprot protein ID}
#'   \item{description}{Protein description}
#'   \item{37C_B1_G1S}{Log2 fold-change G1S vs G1S phase at 37°C in the first bioreplicate}
#'   \item{37C_B1_PM}{Log2 fold-change PM vs G1S phase at 37°C in the first bioreplicate}
#'   \item{37C_B1_S}{Log2 fold-change S vs G1S phase at 37°C in the first bioreplicate}
#'   \item{37C_B2_G1S}{Log2 fold-change G1S vs G1S phase at 37°C in the second bioreplicate}
#'   \item{37C_B2_PM}{Log2 fold-change PM vs G1S phase at 37°C in the second bioreplicate}
#'   \item{37C_B2_S}{Log2 fold-change S vs G1S phase at 37°C in the second bioreplicate}
#'   \item{37C_B3_G1S}{Log2 fold-change G1S vs G1S phase at 37°C in the third bioreplicate}
#'   \item{37C_B3_PM}{Log2 fold-change PM vs G1S phase at 37°C in the third bioreplicate}
#'   \item{37C_B3_S}{Log2 fold-change S vs G1S phase at 37°C in the third bioreplicate}
#'   \item{47C_B1_G1S}{Log2 fold-change G1S vs G1S phase at 47°C in the first bioreplicate}
#'   \item{47C_B1_PM}{Log2 fold-change PM vs G1S phase at 47°C in the first bioreplicate}
#'   \item{47C_B1_S}{Log2 fold-change S vs G1S phase at 47°C in the first bioreplicate}
#'   \item{47C_B2_G1S}{Log2 fold-change G1S vs G1S phase at 47°C in the second bioreplicate}
#'   \item{47C_B2_PM}{Log2 fold-change PM vs G1S phase at 47°C in the second bioreplicate}
#'   \item{47C_B2_S}{Log2 fold-change S vs G1S phase at 47°C in the second bioreplicate}
#'   \item{47C_B3_G1S}{Log2 fold-change G1S vs G1S phase at 47°C in the third bioreplicate}
#'   \item{47C_B3_PM}{Log2 fold-change PM vs G1S phase at 47°C in the third bioreplicate}
#'   \item{47C_B3_S}{Log2 fold-change S vs G1S phase at 47°C in the third bioreplicate}
#'   \item{50C_B1_G1S}{Log2 fold-change G1S vs G1S phase at 50°C in the first bioreplicate}
#'   \item{50C_B1_PM}{Log2 fold-change PM vs G1S phase at 50°C in the first bioreplicate}
#'   \item{50C_B1_S}{Log2 fold-change S vs G1S phase at 50°C in the first bioreplicate}
#'   \item{50C_B2_G1S}{Log2 fold-change G1S vs G1S phase at 50°C in the second bioreplicate}
#'   \item{50C_B2_PM}{Log2 fold-change PM vs G1S phase at 50°C in the second bioreplicate}
#'   \item{50C_B2_S}{Log2 fold-change S vs G1S phase at 50°C in the second bioreplicate}
#'   \item{50C_B3_G1S}{Log2 fold-change G1S vs G1S phase at 50°C in the third bioreplicate}
#'   \item{50C_B3_PM}{Log2 fold-change PM vs G1S phase at 50°C in the third bioreplicate}
#'   \item{50C_B3_S}{Log2 fold-change S vs G1S phase at 50°C in the third bioreplicate}
#'   \item{52C_B1_G1S}{Log2 fold-change G1S vs G1S phase at 52°C in the first bioreplicate}
#'   \item{52C_B1_PM}{Log2 fold-change PM vs G1S phase at 52°C in the first bioreplicate}
#'   \item{52C_B1_S}{Log2 fold-change S vs G1S phase at 52°C in the first bioreplicate}
#'   \item{52C_B2_G1S}{Log2 fold-change G1S vs G1S phase at 52°C in the second bioreplicate}
#'   \item{52C_B2_PM}{Log2 fold-change PM vs G1S phase at 52°C in the second bioreplicate}
#'   \item{52C_B2_S}{Log2 fold-change S vs G1S phase at 52°C in the second bioreplicate}
#'   \item{52C_B3_G1S}{Log2 fold-change G1S vs G1S phase at 52°C in the third bioreplicate}
#'   \item{52C_B3_PM}{Log2 fold-change PM vs G1S phase at 52°C in the third bioreplicate}
#'   \item{52C_B3_S}{Log2 fold-change S vs G1S phase at 52°C in the third bioreplicate}
#'   \item{54C_B1_G1S}{Log2 fold-change G1S vs G1S phase at 54°C in the first bioreplicate}
#'   \item{54C_B1_PM}{Log2 fold-change PM vs G1S phase at 54°C in the first bioreplicate}
#'   \item{54C_B1_S}{Log2 fold-change S vs G1S phase at 54°C in the first bioreplicate}
#'   \item{54C_B2_G1S}{Log2 fold-change G1S vs G1S phase at 54°C in the second bioreplicate}
#'   \item{54C_B2_PM}{Log2 fold-change PM vs G1S phase at 54°C in the second bioreplicate}
#'   \item{54C_B2_S}{Log2 fold-change S vs G1S phase at 54°C in the second bioreplicate}
#'   \item{54C_B3_G1S}{Log2 fold-change G1S vs G1S phase at 54°C in the third bioreplicate}
#'   \item{54C_B3_PM}{Log2 fold-change PM vs G1S phase at 54°C in the third bioreplicate}
#'   \item{54C_B3_S}{Log2 fold-change S vs G1S phase at 54°C in the third bioreplicate}
#'   \item{57C_B1_G1S}{Log2 fold-change G1S vs G1S phase at 57°C in the first bioreplicate}
#'   \item{57C_B1_PM}{Log2 fold-change PM vs G1S phase at 57°C in the first bioreplicate}
#'   \item{57C_B1_S}{Log2 fold-change S vs G1S phase at 57°C in the first bioreplicate}
#'   \item{57C_B2_G1S}{Log2 fold-change G1S vs G1S phase at 57°C in the second bioreplicate}
#'   \item{57C_B2_PM}{Log2 fold-change PM vs G1S phase at 57°C in the second bioreplicate}
#'   \item{57C_B2_S}{Log2 fold-change S vs G1S phase at 57°C in the second bioreplicate}
#'   \item{57C_B3_G1S}{Log2 fold-change G1S vs G1S phase at 57°C in the third bioreplicate}
#'   \item{57C_B3_PM}{Log2 fold-change PM vs G1S phase at 57°C in the third bioreplicate}
#'   \item{57C_B3_S}{Log2 fold-change S vs G1S phase at 57°C in the third bioreplicate}
#'   \item{sumUniPeps}{Median of the number of unique peptides from all temperatures}
#'   \item{sumPSMs}{Median of the PSMs from all temperatures}
#'   \item{countNum}{Median of the abundance count from all reference channel}
#' }
#' @source This data were obtained from 7.	Dai, L., T. Zhao, X. Bisteau et al. (2018)
#'         Modulation of Protein-Interaction States through the Cell Cycle.
#'         Cell, 173, 1481-1494.e13.
"chemarrest"


#' Chemarrest average fold-change dataset
#'
#' An IMPRINTS-CETSA dataset containing the mean log2 fold-change between different cell cycle phase.
#'
#' @format ## `chemarrest_ave`
#' A data frame with 2777 rows and 23 columns:
#' \describe{
#'   \item{id}{Uniprot protein ID}
#'   \item{description}{Protein description}
#'   \item{37C_G1S}{Log2 fold-change G1S vs G1S phase at 37°C}
#'   \item{37C_PM}{Log2 fold-change PM vs G1S phase at 37°C}
#'   \item{37C_S}{Log2 fold-change S vs G1S phase at 37°C}
#'   \item{47C_G1S}{Log2 fold-change G1S vs G1S phase at 47°C}
#'   \item{47C_PM}{Log2 fold-change PM vs G1S phase at 47°C}
#'   \item{47C_S}{Log2 fold-change S vs G1S phase at 47°C}
#'   \item{50C_G1S}{Log2 fold-change G1S vs G1S phase at 50°C}
#'   \item{50C_PM}{Log2 fold-change PM vs G1S phase at 50°C}
#'   \item{50C_S}{Log2 fold-change S vs G1S phase at 50°C}
#'   \item{52C_G1S}{Log2 fold-change G1S vs G1S phase at 52°C}
#'   \item{52C_PM}{Log2 fold-change PM vs G1S phase at 52°C}
#'   \item{52C_S}{Log2 fold-change S vs G1S phase at 52°C}
#'   \item{54C_G1S}{Log2 fold-change G1S vs G1S phase at 54°C}
#'   \item{54C_PM}{Log2 fold-change PM vs G1S phase at 54°C}
#'   \item{54C_S}{Log2 fold-change S vs G1S phase at 54°C}
#'   \item{57C_G1S}{Log2 fold-change G1S vs G1S phase at 57°C}
#'   \item{57C_PM}{Log2 fold-change PM vs G1S phase at 57°C}
#'   \item{57C_S}{Log2 fold-change S vs G1S phase at 57°C}
#'   \item{sumUniPeps}{Median of the number of unique peptides from all temperatures}
#'   \item{sumPSMs}{Median of the PSMs from all temperatures}
#'   \item{countNum}{Median of the abundance count from all reference channel}
#' }
#' @source This data were obtained from 7.	Dai, L., T. Zhao, X. Bisteau et al. (2018)
#'         Modulation of Protein-Interaction States through the Cell Cycle.
#'         Cell, 173, 1481-1494.e13.
"chemarrest_ave"


#' Hitlist from the chemarrest dataset
#'
#' An hitlist from the chemarrest dataset obtained with the \code{imprints_score} function
#' from the package \code{IMPRINTS.CETSA} using default parameters except useMAD set to TRUE.
#'
#' @format ## `hitlist_chemarrest`
#' A data frame with 732 rows and 3 columns:
#' \describe{
#'   \item{id}{Uniprot protein ID}
#'   \item{treatment}{The treatment in which it was identified as a hit}
#'   \item{category}{Category of the hit}
#' }
#' @source This data were obtained from 7.	Dai, L., T. Zhao, X. Bisteau et al. (2018)
#'         Modulation of Protein-Interaction States through the Cell Cycle.
#'         Cell, 173, 1481-1494.e13.
"hitlist_chemarrest"


#' The non hits from the chemarrest dataset
#'
#' The non hits derived from the hitlist of the chemarrest dataset
#'
#' @format ## `NN_chemarrest`
#' A data frame with 4706 rows and 4 columns:
#' \describe{
#'   \item{id}{Uniprot protein ID}
#'   \item{description}{Protein description}
#'   \item{treatment}{The treatment in which it was identified as a hit}
#'   \item{category}{Category of the non hit, i.e. NN}
#' }
#' @source This data were obtained from 7.	Dai, L., T. Zhao, X. Bisteau et al. (2018)
#'         Modulation of Protein-Interaction States through the Cell Cycle.
#'         Cell, 173, 1481-1494.e13.
"NN_chemarrest"


#' Raw Protein.Group from the shortened  elutriation dataset
#'
#' A data frame containing the six Protein.Group files from the raw data of the elutriation dataset.
#' These raw data have been processed with Proteome Discovers version 2.4 and shortened  to
#' approximately 500 proteins to showcase the functionnalities of the app and show the expected format.
#'
#' @format ## `elu_example1_raw_joined`
#' A data frame with 3170 rows and 16 columns:
#' \describe{
#'   \item{id}{Uniprot protein ID}
#'   \item{description}{Protein description}
#'   \item{condition}{}
#'   \item{B1_G1}{Raw protein abundance of the first bioreplicate of the G1 phase}
#'   \item{B1_S}{Raw protein abundance of the first bioreplicate of the S phase}
#'   \item{B1_G2}{Raw protein abundance of the first bioreplicate of the G2 phase}
#'   \item{B2_G1}{Raw protein abundance of the second bioreplicate of the G1 phase}
#'   \item{B2_S}{Raw protein abundance of the second bioreplicate of the S phase}
#'   \item{B2_G2}{Raw protein abundance of the second bioreplicate of the G2 phase}
#'   \item{B3_S}{Raw protein abundance of the third bioreplicate of the S phase}
#'   \item{B3_G1}{Raw protein abundance of the third bioreplicate of the G1 phase}
#'   \item{B3_G2}{Raw protein abundance of the third bioreplicate of the G2 phase}
#'   \item{Mix}{Raw protein abundance of the Mix channel}
#'   \item{sumUniPeps}{Number of unique peptides}
#'   \item{sumPSMs}{PSMs}
#'   \item{countNum}{Abundance count from reference channel}
#' }
#' @source This data were obtained from 7.	Dai, L., T. Zhao, X. Bisteau et al. (2018)
#'         Modulation of Protein-Interaction States through the Cell Cycle.
#'         Cell, 173, 1481-1494.e13.
"elu_example1_raw_joined"


#' Filtered Protein.Group from the shortened  elutriation dataset
#'
#' A data frame containing the cleaned and filtered raw protein abundance using the \code{ms_clean} function
#' from \code{IMPRINTS.CETSA} package.
#'
#' See \code{\link{elu_example1_raw_joined}}
#'
#' @format ## `elu_example2_raw_cleaned`
#' A data frame with 3167 rows and 15 columns:
#' \describe{
#'   \item{id}{Uniprot protein ID}
#'   \item{description}{Protein description}
#'   \item{condition}{}
#'   \item{B1_G1}{Raw protein abundance of the first bioreplicate of the G1 phase}
#'   \item{B1_S}{Raw protein abundance of the first bioreplicate of the S phase}
#'   \item{B1_G2}{Raw protein abundance of the first bioreplicate of the G2 phase}
#'   \item{B2_G1}{Raw protein abundance of the second bioreplicate of the G1 phase}
#'   \item{B2_S}{Raw protein abundance of the second bioreplicate of the S phase}
#'   \item{B2_G2}{Raw protein abundance of the second bioreplicate of the G2 phase}
#'   \item{B3_S}{Raw protein abundance of the third bioreplicate of the S phase}
#'   \item{B3_G1}{Raw protein abundance of the third bioreplicate of the G1 phase}
#'   \item{B3_G2}{Raw protein abundance of the third bioreplicate of the G2 phase}
#'   \item{sumUniPeps}{Number of unique peptides}
#'   \item{sumPSMs}{PSMs}
#'   \item{countNum}{Abundance count from reference channel}
#' }
#' @source This data were obtained from 7.	Dai, L., T. Zhao, X. Bisteau et al. (2018)
#'         Modulation of Protein-Interaction States through the Cell Cycle.
#'         Cell, 173, 1481-1494.e13.
"elu_example2_raw_cleaned"

#' Isorform resolved Protein.Group from the shortened  elutriation dataset
#'
#' The cleaned and filtered raw protein abundance with solved isoforms
#' using  \code{ms_isoform_resolve} function from \code{IMPRINTS.CETSA} package.
#'
#' See \code{\link{elu_example1_raw_joined}}
#'
#' @format ## `elu_example3_isoform_resolved`
#' A data frame with 3167 rows and 15 columns:
#' \describe{
#'   \item{id}{Uniprot protein ID}
#'   \item{description}{Protein description}
#'   \item{condition}{}
#'   \item{B1_G1}{Raw protein abundance of the first bioreplicate of the G1 phase}
#'   \item{B1_S}{Raw protein abundance of the first bioreplicate of the S phase}
#'   \item{B1_G2}{Raw protein abundance of the first bioreplicate of the G2 phase}
#'   \item{B2_G1}{Raw protein abundance of the second bioreplicate of the G1 phase}
#'   \item{B2_S}{Raw protein abundance of the second bioreplicate of the S phase}
#'   \item{B2_G2}{Raw protein abundance of the second bioreplicate of the G2 phase}
#'   \item{B3_S}{Raw protein abundance of the third bioreplicate of the S phase}
#'   \item{B3_G1}{Raw protein abundance of the third bioreplicate of the G1 phase}
#'   \item{B3_G2}{Raw protein abundance of the third bioreplicate of the G2 phase}
#'   \item{sumUniPeps}{Number of unique peptides}
#'   \item{sumPSMs}{PSMs}
#'   \item{countNum}{Abundance count from reference channel}
#' }
#' @source This data were obtained from 7.	Dai, L., T. Zhao, X. Bisteau et al. (2018)
#'         Modulation of Protein-Interaction States through the Cell Cycle.
#'         Cell, 173, 1481-1494.e13.
"elu_example3_isoform_resolved"

#' Rearranged protein abundance from the shortened elutriation dataset
#'
#' A data frame containing the raw abundance protein from the shortened elutriation dataset put in the right
#' IMPRINTS-CETSA format.
#'
#' See \code{\link{elu_example1_raw_joined}}
#'
#' @format ## `elu_example4_pre_norm`
#' A data frame with 521 rows and 59 columns:
#' \describe{
#'   \item{id}{Uniprot protein ID}
#'   \item{description}{Protein description}
#'   \item{37C_B1_G1}{Raw protein abundance in G1 phase at 37°C in the first bioreplicate}
#'   \item{37C_B1_G2}{Raw protein abundance in G2 phase at 37°C in the first bioreplicate}
#'   \item{37C_B1_S}{Raw protein abundance in S phase at 37°C in the first bioreplicate}
#'   \item{37C_B2_G1}{Raw protein abundance in G1 phase at 37°C in the second bioreplicate}
#'   \item{37C_B2_G2}{Raw protein abundance in G2 phase at 37°C in the second bioreplicate}
#'   \item{37C_B2_S}{Raw protein abundance in S phase at 37°C in the second bioreplicate}
#'   \item{37C_B3_G1}{Raw protein abundance in G1 phase at 37°C in the third bioreplicate}
#'   \item{37C_B3_G2}{Raw protein abundance in G2 phase at 37°C in the third bioreplicate}
#'   \item{37C_B3_S}{Raw protein abundance in S phase at 37°C in the third bioreplicate}
#'   \item{47C_B1_G1}{Raw protein abundance in G1 phase at 47°C in the first bioreplicate}
#'   \item{47C_B1_G2}{Raw protein abundance in G2 phase at 47°C in the first bioreplicate}
#'   \item{47C_B1_S}{Raw protein abundance in S phase at 47°C in the first bioreplicate}
#'   \item{47C_B2_G1}{Raw protein abundance in G1 phase at 47°C in the second bioreplicate}
#'   \item{47C_B2_G2}{Raw protein abundance in G2 phase at 47°C in the second bioreplicate}
#'   \item{47C_B2_S}{Raw protein abundance in S phase at 47°C in the second bioreplicate}
#'   \item{47C_B3_G1}{Raw protein abundance in G1 phase at 47°C in the third bioreplicate}
#'   \item{47C_B3_G2}{Raw protein abundance in G2 phase at 47°C in the third bioreplicate}
#'   \item{47C_B3_S}{Raw protein abundance in S phase at 47°C in the third bioreplicate}
#'   \item{50C_B1_G1}{Raw protein abundance in G1 phase at 50°C in the first bioreplicate}
#'   \item{50C_B1_G2}{Raw protein abundance in G2 phase at 50°C in the first bioreplicate}
#'   \item{50C_B1_S}{Raw protein abundance in S phase at 50°C in the first bioreplicate}
#'   \item{50C_B2_G1}{Raw protein abundance in G1 phase at 50°C in the second bioreplicate}
#'   \item{50C_B2_G2}{Raw protein abundance in G2 phase at 50°C in the second bioreplicate}
#'   \item{50C_B2_S}{Raw protein abundance in S phase at 50°C in the second bioreplicate}
#'   \item{50C_B3_G1}{Raw protein abundance in G1 phase at 50°C in the third bioreplicate}
#'   \item{50C_B3_G2}{Raw protein abundance in G2 phase at 50°C in the third bioreplicate}
#'   \item{50C_B3_S}{Raw protein abundance in S phase at 50°C in the third bioreplicate}
#'   \item{52C_B1_G1}{Raw protein abundance in G1 phase at 52°C in the first bioreplicate}
#'   \item{52C_B1_G2}{Raw protein abundance in G2 phase at 52°C in the first bioreplicate}
#'   \item{52C_B1_S}{Raw protein abundance in S phase at 52°C in the first bioreplicate}
#'   \item{52C_B2_G1}{Raw protein abundance in G1 phase at 52°C in the second bioreplicate}
#'   \item{52C_B2_G2}{Raw protein abundance in G2 phase at 52°C in the second bioreplicate}
#'   \item{52C_B2_S}{Raw protein abundance in S phase at 52°C in the second bioreplicate}
#'   \item{52C_B3_G1}{Raw protein abundance in G1 phase at 52°C in the third bioreplicate}
#'   \item{52C_B3_G2}{Raw protein abundance in G2 phase at 52°C in the third bioreplicate}
#'   \item{52C_B3_S}{Raw protein abundance in S phase at 52°C in the third bioreplicate}
#'   \item{54C_B1_G1}{Raw protein abundance in G1 phase at 54°C in the first bioreplicate}
#'   \item{54C_B1_G2}{Raw protein abundance in G2 phase at 54°C in the first bioreplicate}
#'   \item{54C_B1_S}{Raw protein abundance in S phase at 54°C in the first bioreplicate}
#'   \item{54C_B2_G1}{Raw protein abundance in G1 phase at 54°C in the second bioreplicate}
#'   \item{54C_B2_G2}{Raw protein abundance in G2 phase at 54°C in the second bioreplicate}
#'   \item{54C_B2_S}{Raw protein abundance in S phase at 54°C in the second bioreplicate}
#'   \item{54C_B3_G1}{Raw protein abundance in G1 phase at 54°C in the third bioreplicate}
#'   \item{54C_B3_G2}{Raw protein abundance in G2 phase at 54°C in the third bioreplicate}
#'   \item{54C_B3_S}{Raw protein abundance in S phase at 54°C in the third bioreplicate}
#'   \item{57C_B1_G1}{Raw protein abundance in G1 phase at 57°C in the first bioreplicate}
#'   \item{57C_B1_G2}{Raw protein abundance in G2 phase at 57°C in the first bioreplicate}
#'   \item{57C_B1_S}{Raw protein abundance in S phase at 57°C in the first bioreplicate}
#'   \item{57C_B2_G1}{Raw protein abundance in G1 phase at 57°C in the second bioreplicate}
#'   \item{57C_B2_G2}{Raw protein abundance in G2 phase at 57°C in the second bioreplicate}
#'   \item{57C_B2_S}{Raw protein abundance in S phase at 57°C in the second bioreplicate}
#'   \item{57C_B3_G1}{Raw protein abundance in G1 phase at 57°C in the third bioreplicate}
#'   \item{57C_B3_G2}{Raw protein abundance in G2 phase at 57°C in the third bioreplicate}
#'   \item{57C_B3_S}{Raw protein abundance in S phase at 57°C in the third bioreplicate}
#'   \item{sumUniPeps}{Median of the number of unique peptides from all temperatures}
#'   \item{sumPSMs}{Median of the PSMs from all temperatures}
#'   \item{countNum}{Median of the abundance count from all reference channel}
#' }
#' @source This data were obtained from 7.	Dai, L., T. Zhao, X. Bisteau et al. (2018)
#'         Modulation of Protein-Interaction States through the Cell Cycle.
#'         Cell, 173, 1481-1494.e13.
"elu_example4_pre_norm"

#' Normalized protein abundance from the shortened  elutriation dataset
#'
#' A data frame containing the normalized abundance protein from the shortened
#' elutriation dataset.
#'
#' See \code{\link{elu_example1_raw_joined}}
#'
#' @format ## `elu_example5_post_norm`
#' A data frame with 521 rows and 59 columns:
#' \describe{
#'   \item{id}{Uniprot protein ID}
#'   \item{description}{Protein description}
#'   \item{37C_B1_G1}{Normalized protein abundance in G1 phase at 37°C in the first bioreplicate}
#'   \item{37C_B1_G2}{Normalized protein abundance in G2 phase at 37°C in the first bioreplicate}
#'   \item{37C_B1_S}{Normalized protein abundance in S phase at 37°C in the first bioreplicate}
#'   \item{37C_B2_G1}{Normalized protein abundance in G1 phase at 37°C in the second bioreplicate}
#'   \item{37C_B2_G2}{Normalized protein abundance in G2 phase at 37°C in the second bioreplicate}
#'   \item{37C_B2_S}{Normalized protein abundance in S phase at 37°C in the second bioreplicate}
#'   \item{37C_B3_G1}{Normalized protein abundance in G1 phase at 37°C in the third bioreplicate}
#'   \item{37C_B3_G2}{Normalized protein abundance in G2 phase at 37°C in the third bioreplicate}
#'   \item{37C_B3_S}{Normalized protein abundance in S phase at 37°C in the third bioreplicate}
#'   \item{47C_B1_G1}{Normalized protein abundance in G1 phase at 47°C in the first bioreplicate}
#'   \item{47C_B1_G2}{Normalized protein abundance in G2 phase at 47°C in the first bioreplicate}
#'   \item{47C_B1_S}{Normalized protein abundance in S phase at 47°C in the first bioreplicate}
#'   \item{47C_B2_G1}{Normalized protein abundance in G1 phase at 47°C in the second bioreplicate}
#'   \item{47C_B2_G2}{Normalized protein abundance in G2 phase at 47°C in the second bioreplicate}
#'   \item{47C_B2_S}{Normalized protein abundance in S phase at 47°C in the second bioreplicate}
#'   \item{47C_B3_G1}{Normalized protein abundance in G1 phase at 47°C in the third bioreplicate}
#'   \item{47C_B3_G2}{Normalized protein abundance in G2 phase at 47°C in the third bioreplicate}
#'   \item{47C_B3_S}{Normalized protein abundance in S phase at 47°C in the third bioreplicate}
#'   \item{50C_B1_G1}{Normalized protein abundance in G1 phase at 50°C in the first bioreplicate}
#'   \item{50C_B1_G2}{Normalized protein abundance in G2 phase at 50°C in the first bioreplicate}
#'   \item{50C_B1_S}{Normalized protein abundance in S phase at 50°C in the first bioreplicate}
#'   \item{50C_B2_G1}{Normalized protein abundance in G1 phase at 50°C in the second bioreplicate}
#'   \item{50C_B2_G2}{Normalized protein abundance in G2 phase at 50°C in the second bioreplicate}
#'   \item{50C_B2_S}{Normalized protein abundance in S phase at 50°C in the second bioreplicate}
#'   \item{50C_B3_G1}{Normalized protein abundance in G1 phase at 50°C in the third bioreplicate}
#'   \item{50C_B3_G2}{Normalized protein abundance in G2 phase at 50°C in the third bioreplicate}
#'   \item{50C_B3_S}{Normalized protein abundance in S phase at 50°C in the third bioreplicate}
#'   \item{52C_B1_G1}{Normalized protein abundance in G1 phase at 52°C in the first bioreplicate}
#'   \item{52C_B1_G2}{Normalized protein abundance in G2 phase at 52°C in the first bioreplicate}
#'   \item{52C_B1_S}{Normalized protein abundance in S phase at 52°C in the first bioreplicate}
#'   \item{52C_B2_G1}{Normalized protein abundance in G1 phase at 52°C in the second bioreplicate}
#'   \item{52C_B2_G2}{Normalized protein abundance in G2 phase at 52°C in the second bioreplicate}
#'   \item{52C_B2_S}{Normalized protein abundance in S phase at 52°C in the second bioreplicate}
#'   \item{52C_B3_G1}{Normalized protein abundance in G1 phase at 52°C in the third bioreplicate}
#'   \item{52C_B3_G2}{Normalized protein abundance in G2 phase at 52°C in the third bioreplicate}
#'   \item{52C_B3_S}{Normalized protein abundance in S phase at 52°C in the third bioreplicate}
#'   \item{54C_B1_G1}{Normalized protein abundance in G1 phase at 54°C in the first bioreplicate}
#'   \item{54C_B1_G2}{Normalized protein abundance in G2 phase at 54°C in the first bioreplicate}
#'   \item{54C_B1_S}{Normalized protein abundance in S phase at 54°C in the first bioreplicate}
#'   \item{54C_B2_G1}{Normalized protein abundance in G1 phase at 54°C in the second bioreplicate}
#'   \item{54C_B2_G2}{Normalized protein abundance in G2 phase at 54°C in the second bioreplicate}
#'   \item{54C_B2_S}{Normalized protein abundance in S phase at 54°C in the second bioreplicate}
#'   \item{54C_B3_G1}{Normalized protein abundance in G1 phase at 54°C in the third bioreplicate}
#'   \item{54C_B3_G2}{Normalized protein abundance in G2 phase at 54°C in the third bioreplicate}
#'   \item{54C_B3_S}{Normalized protein abundance in S phase at 54°C in the third bioreplicate}
#'   \item{57C_B1_G1}{Normalized protein abundance in G1 phase at 57°C in the first bioreplicate}
#'   \item{57C_B1_G2}{Normalized protein abundance in G2 phase at 57°C in the first bioreplicate}
#'   \item{57C_B1_S}{Normalized protein abundance in S phase at 57°C in the first bioreplicate}
#'   \item{57C_B2_G1}{Normalized protein abundance in G1 phase at 57°C in the second bioreplicate}
#'   \item{57C_B2_G2}{Normalized protein abundance in G2 phase at 57°C in the second bioreplicate}
#'   \item{57C_B2_S}{Normalized protein abundance in S phase at 57°C in the second bioreplicate}
#'   \item{57C_B3_G1}{Normalized protein abundance in G1 phase at 57°C in the third bioreplicate}
#'   \item{57C_B3_G2}{Normalized protein abundance in G2 phase at 57°C in the third bioreplicate}
#'   \item{57C_B3_S}{Normalized protein abundance in S phase at 57°C in the third bioreplicate}
#'   \item{sumUniPeps}{Median of the number of unique peptides from all temperatures}
#'   \item{sumPSMs}{Median of the PSMs from all temperatures}
#'   \item{countNum}{Median of the abundance count from all reference channel}
#' }
#' @source This data were obtained from 7.	Dai, L., T. Zhao, X. Bisteau et al. (2018)
#'         Modulation of Protein-Interaction States through the Cell Cycle.
#'         Cell, 173, 1481-1494.e13.
"elu_example5_post_norm"


#' Log2 fold change from the shortened  elutriation dataset
#'
#' A data frame containing the log2 fold change from the shortened
#' elutriation dataset.
#'
#' See \code{\link{elu_example1_raw_joined}}
#'
#' @format ## `elu_example6_caldiff`
#' A data frame with 521 rows and 59 columns:
#' \describe{
#'   \item{id}{Uniprot protein ID}
#'   \item{description}{Protein description}
#'   \item{37C_B1_G1}{Log2 fold-change G1 vs G1 phase at 37°C in the first bioreplicate}
#'   \item{37C_B1_G2}{Log2 fold-change G2 vs G1 phase at 37°C in the first bioreplicate}
#'   \item{37C_B1_S}{Log2 fold-change S vs G1 phase at 37°C in the first bioreplicate}
#'   \item{37C_B2_G1}{Log2 fold-change G1 vs G1 phase at 37°C in the second bioreplicate}
#'   \item{37C_B2_G2}{Log2 fold-change G2 vs G1 phase at 37°C in the second bioreplicate}
#'   \item{37C_B2_S}{Log2 fold-change S vs G1 phase at 37°C in the second bioreplicate}
#'   \item{37C_B3_G1}{Log2 fold-change G1 vs G1 phase at 37°C in the third bioreplicate}
#'   \item{37C_B3_G2}{Log2 fold-change G2 vs G1 phase at 37°C in the third bioreplicate}
#'   \item{37C_B3_S}{Log2 fold-change S vs G1 phase at 37°C in the third bioreplicate}
#'   \item{47C_B1_G1}{Log2 fold-change G1 vs G1 phase at 47°C in the first bioreplicate}
#'   \item{47C_B1_G2}{Log2 fold-change G2 vs G1 phase at 47°C in the first bioreplicate}
#'   \item{47C_B1_S}{Log2 fold-change S vs G1 phase at 47°C in the first bioreplicate}
#'   \item{47C_B2_G1}{Log2 fold-change G1 vs G1 phase at 47°C in the second bioreplicate}
#'   \item{47C_B2_G2}{Log2 fold-change G2 vs G1 phase at 47°C in the second bioreplicate}
#'   \item{47C_B2_S}{Log2 fold-change S vs G1 phase at 47°C in the second bioreplicate}
#'   \item{47C_B3_G1}{Log2 fold-change G1 vs G1 phase at 47°C in the third bioreplicate}
#'   \item{47C_B3_G2}{Log2 fold-change G2 vs G1 phase at 47°C in the third bioreplicate}
#'   \item{47C_B3_S}{Log2 fold-change S vs G1 phase at 47°C in the third bioreplicate}
#'   \item{50C_B1_G1}{Log2 fold-change G1 vs G1 phase at 50°C in the first bioreplicate}
#'   \item{50C_B1_G2}{Log2 fold-change G2 vs G1 phase at 50°C in the first bioreplicate}
#'   \item{50C_B1_S}{Log2 fold-change S vs G1 phase at 50°C in the first bioreplicate}
#'   \item{50C_B2_G1}{Log2 fold-change G1 vs G1 phase at 50°C in the second bioreplicate}
#'   \item{50C_B2_G2}{Log2 fold-change G2 vs G1 phase at 50°C in the second bioreplicate}
#'   \item{50C_B2_S}{Log2 fold-change S vs G1 phase at 50°C in the second bioreplicate}
#'   \item{50C_B3_G1}{Log2 fold-change G1 vs G1 phase at 50°C in the third bioreplicate}
#'   \item{50C_B3_G2}{Log2 fold-change G2 vs G1 phase at 50°C in the third bioreplicate}
#'   \item{50C_B3_S}{Log2 fold-change S vs G1 phase at 50°C in the third bioreplicate}
#'   \item{52C_B1_G1}{Log2 fold-change G1 vs G1 phase at 52°C in the first bioreplicate}
#'   \item{52C_B1_G2}{Log2 fold-change G2 vs G1 phase at 52°C in the first bioreplicate}
#'   \item{52C_B1_S}{Log2 fold-change S vs G1 phase at 52°C in the first bioreplicate}
#'   \item{52C_B2_G1}{Log2 fold-change G1 vs G1 phase at 52°C in the second bioreplicate}
#'   \item{52C_B2_G2}{Log2 fold-change G2 vs G1 phase at 52°C in the second bioreplicate}
#'   \item{52C_B2_S}{Log2 fold-change S vs G1 phase at 52°C in the second bioreplicate}
#'   \item{52C_B3_G1}{Log2 fold-change G1 vs G1 phase at 52°C in the third bioreplicate}
#'   \item{52C_B3_G2}{Log2 fold-change G2 vs G1 phase at 52°C in the third bioreplicate}
#'   \item{52C_B3_S}{Log2 fold-change S vs G1 phase at 52°C in the third bioreplicate}
#'   \item{54C_B1_G1}{Log2 fold-change G1 vs G1 phase at 54°C in the first bioreplicate}
#'   \item{54C_B1_G2}{Log2 fold-change G2 vs G1 phase at 54°C in the first bioreplicate}
#'   \item{54C_B1_S}{Log2 fold-change S vs G1 phase at 54°C in the first bioreplicate}
#'   \item{54C_B2_G1}{Log2 fold-change G1 vs G1 phase at 54°C in the second bioreplicate}
#'   \item{54C_B2_G2}{Log2 fold-change G2 vs G1 phase at 54°C in the second bioreplicate}
#'   \item{54C_B2_S}{Log2 fold-change S vs G1 phase at 54°C in the second bioreplicate}
#'   \item{54C_B3_G1}{Log2 fold-change G1 vs G1 phase at 54°C in the third bioreplicate}
#'   \item{54C_B3_G2}{Log2 fold-change G2 vs G1 phase at 54°C in the third bioreplicate}
#'   \item{54C_B3_S}{Log2 fold-change S vs G1 phase at 54°C in the third bioreplicate}
#'   \item{57C_B1_G1}{Log2 fold-change G1 vs G1 phase at 57°C in the first bioreplicate}
#'   \item{57C_B1_G2}{Log2 fold-change G2 vs G1 phase at 57°C in the first bioreplicate}
#'   \item{57C_B1_S}{Log2 fold-change S vs G1 phase at 57°C in the first bioreplicate}
#'   \item{57C_B2_G1}{Log2 fold-change G1 vs G1 phase at 57°C in the second bioreplicate}
#'   \item{57C_B2_G2}{Log2 fold-change G2 vs G1 phase at 57°C in the second bioreplicate}
#'   \item{57C_B2_S}{Log2 fold-change S vs G1 phase at 57°C in the second bioreplicate}
#'   \item{57C_B3_G1}{Log2 fold-change G1 vs G1 phase at 57°C in the third bioreplicate}
#'   \item{57C_B3_G2}{Log2 fold-change G2 vs G1 phase at 57°C in the third bioreplicate}
#'   \item{57C_B3_S}{Log2 fold-change S vs G1 phase at 57°C in the third bioreplicate}
#'   \item{sumUniPeps}{Median of the number of unique peptides from all temperatures}
#'   \item{sumPSMs}{Median of the PSMs from all temperatures}
#'   \item{countNum}{Median of the abundance count from all reference channel}
#' }
#' @source This data were obtained from 7.	Dai, L., T. Zhao, X. Bisteau et al. (2018)
#'         Modulation of Protein-Interaction States through the Cell Cycle.
#'         Cell, 173, 1481-1494.e13.
"elu_example6_caldiff"


#' Human protein cellular localization from Protein Atlas database
#'
#'
#' @format ## `pr_atlas`
#' A data frame with 17 377 rows and 6 columns:
#' \describe{
#'   \item{Gene}{Gene}
#'   \item{Uniprot}{Uniprot ID}
#'   \item{Reliability..IH.}{Reliability..IH.}
#'   \item{Subcellular.location}{Subcellular.location}
#'   \item{Subcellular.main.location}{Subcellular.main.location}
#'   \item{Subcellular.additional.location}{Subcellular.additional.location}
#' }
#' @source https://www.proteinatlas.org/api/search_download.php?search=Human&format=json&columns=g,up,relih,scl,scml,scal&compress=no
#'   Last downloaded on the 21th of December 2023
"pr_atlas"


#' Mouse protein cellular localization from Protein Atlas database
#'
#'
#' @format ## `pr_atlas_mouse`
#' A data frame with 1525 rows and 6 columns:
#' \describe{
#'   \item{Gene}{Gene}
#'   \item{Uniprot}{Uniprot ID}
#'   \item{Reliability..IH.}{Reliability..IH.}
#'   \item{Subcellular.location}{Subcellular.location}
#'   \item{Subcellular.main.location}{Subcellular.main.location}
#'   \item{Subcellular.additional.location}{Subcellular.additional.location}
#' }
#' @source https://www.proteinatlas.org/api/search_download.php?search=Mouse&format=json&columns=g,up,relih,scl,scml,scal&compress=no
#'   Last downloaded on the 21th of December 2023
"pr_atlas_mouse"
