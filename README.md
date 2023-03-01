# mineCETSAapp
This is a package to complete the IMRPRINTS.CETSA package from Dai Lingyun (https://github.com/nkdailingyun/IMRPRINTS.CETSA), more specifically the CETSA IMPRINTS analysis.
The goal of this package is to provide a user friendly interface of this package. With this shiny app, you can do the complete IMPRINTS CETSA analysis on the protein and peptide level; visualize the 2D bar plot or heatmap of many proteins from drug experiment available in the database; add new dataset or remove old ones from this database; get STRING network and extract other GO analysis from your data; plot interactive barplots networks; start a PubMed search based on the protein names or any character you want and save the results in word files; etc.

The app contains also a tab called "Interactive cell". The goal of this one is to visualize proteins directly into the cell and get information on it when clicking one. The subcellular location is obtained thanks to the protein Atlas data base.

Moreover, the package contains a new way of getting hits with the function called SR_CetsaHit (SR = Stability Rate). This function computes a score and a p-value for each protein and return a volcano plot with the selected FDR. Also, thnaks to this scoring, the proteins can be ranked and used for further gene ontology analysis.
mineCETSAapp also allows to do the IMPRINTS analysis on the peptide level and search for potential cleaved peptides. 

Some data are available in the package so that you can play with when you first try the app. The data are from the [The cell cycle paper (2018)](https://www-sciencedirect-com.proxy.insermbiblio.inist.fr/science/article/pii/S0092867418303970?via%3Dihub) from Dai Lingyun and al.

## What is CETSA?
The Cellular Thermal Shift Assay (CETSA) (orginially described in [Science 341(6141):84-87 (2013)](http://www.sciencemag.org/lookup/doi/10.1126/science.1233606)) is a 
biophysical assay based on the principle of ligand-induced thermal stabilization of target proteins, meaning that a protein's melting temperature will change upon 
ligand interaction.
 
By heating samples (lysate, cells or tissue pieces) to different temperatures, and quantifying proteins in the soluble fraction we can detect altered protein interactions 
after for example drug treatment. This can either be done for selected proteins of interest by using antibody-based assays or on a proteome-wide level by using 
mass spectrometry.  

CETSA allows direct monitoring of ligand binding to a specific target (target engagement) in lysate, live cells or even tissue pieces. 
It can also be used to study downstream effects on protein interaction, providing a novel perspective on protein function in situ. For more details 
please refer to the [CETSA website](https://www.cetsa.org/about). 

This package is specifically conceived to analyse the IMPRINTS-CETSA data. With this technique, less temperatures are required; only fold change bar plot are computed not the melt curve. With this informative technique, you can compare different cell state by temperature in one experiment set.
 

## How to install mineCETSAapp?  
First go to Rstudio. This package totally depends on the last version of the mineCETSA package. If you havn't donwload it yet, please do so.
For installation from github https://github.com/nkdailingyun/mineCETSA ; however it may not be the last version.
Moreover, you need to install some packages from BioConductor especially for the gene ontology analysis :

```c
if(!requireNamespace("BiocManager", quietly = TRUE)){
   install.packages("BiocManager") 
}
BiocManager::install(c("EBImage", "STRINGdb", "clusterProfiler", "biomaRt", "enrichplot", "multtest"))
```

When all of this is done, type and run the following commands in R console:

```c
if(!requireNamespace("devtools", quietly = TRUE)){
   install.packages("devtools")
} 
devtools::install_github("mgerault/mineCETSAapp")

library(mineCETSAapp)
```

If you want now to use the app just type :

```c
runCETSAapp()
```
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
