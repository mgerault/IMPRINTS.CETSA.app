# What is IMPRINTS.CETSA.app ?
This is a package to complete and provide a user friendly interface for the [IMRPRINTS.CETSA package](https://github.com/nkdailingyun/IMPRINTS.CETSA), more specifically for the analysis of the IMPRINTS-CETSA datasets.
With its Shiny app, you can carry out the complete IMPRINTS-CETSA analysis at the protein and peptide level; visualize the bar plot or heatmap of proteins-of-interest; import new datasets or remove old ones at your end; get STRING network and perform GO analysis on your data; plot interactive networks of bar plots; start a PubMed search based on keywords such as protein names and save the results in Word files; etc.

The app contains also a tab called "Interactive cell". The goal of this one is to visualize proteins directly into the cell and get information on it when clicking one. The subcellular location is obtained thanks to the protein Atlas data base.

Moreover, the package contains a new way of getting hits with the function called imprints_IS (IS = Intercept Score). This function computes a score and a p-value for each protein and return a volcano plot with the selected FDR. Also, thanks to this scoring, the proteins can be ranked and used for further gene ontology analysis.
IMPRINTS.CETSA.app also allows to do the IMPRINTS analysis on the peptide level and search for potential cleaved peptides. 

Some Demo data are included in the package for you to play with. The data are from the [The cell cycle paper (2018)](https://www-sciencedirect-com.proxy.insermbiblio.inist.fr/science/article/pii/S0092867418303970?via%3Dihub) from Dai et al.

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
 

## How to install IMPRINTS.CETSA.app ?  
First go to Rstudio. This package totally depends on the last version of the IMRPRINTS.CETSA package. If you havn't donwload it yet, please do so.
For installation from github https://github.com/nkdailingyun/IMPRINTS.CETSA .
Moreover, you need to install some packages from BioConductor especially for the gene ontology analysis :

```c
if(!requireNamespace("BiocManager", quietly = TRUE)){
   install.packages("BiocManager") 
}
BiocManager::install(c("STRINGdb", "clusterProfiler", "biomaRt", "enrichplot", "multtest"))
```

When all of this is done, type and run the following commands in R console:

```c
if(!requireNamespace("devtools", quietly = TRUE)){
   install.packages("devtools")
} 
devtools::install_github("mgerault/IMPRINTS.CETSA.app")

library(IMPRINTS.CETSA.app)
```

If you want now to use the app just type :

```c
runCETSAapp()
```
 
For more information on how to use the app, go check the [tutorial](https://youtu.be/djpP8nc_JUE) video !
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
