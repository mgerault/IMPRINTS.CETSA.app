# mineCETSAapp
This is a package to complete the mineCETSA package from Lingyun nkdailingyun (https://github.com/nkdailingyun/mineCETSA). 
The goal of this package is to provide a user friendly interface of this package. With this shiny app, you can do the complete CETSA 2D analysis; visualize the 2D bar plot
of many proteins from drug experiement available in the database; add new dataset or remove old ones from this database; start a PubMed search based on the protein names or
any character you want and save the results in word files; etc.

The app contains also a tab called "Interactive cell". This one is not perfectly finished but the goal is to visualize proteins directly into the cell
and get information on it when clicking or select one.

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
 

## How to install mineCETSAapp?  
First go to Rstudio. This package totally depends on the last version of the mineCETSA package. If you havn't donwload it yet, please do so.
For installation from github https://github.com/nkdailingyun/mineCETSA ; it may not be the last version also.
Moreover, you need to install EBImage and STRINGdb package from BioConductor :

*> if(!requireNamespace("BiocManager", quietly = TRUE)){*

*> install.packages("BiocManager")*  

*> }*

*> BiocManager::install(c("EBImage", "STRINGdb"))*  

When all of this is done, type and run the following commands in R console:

*> if(!requireNamespace("devtools", quietly = TRUE)){*

*> install.packages("devtools")*  

*> }*  

*> devtools::install_github("mgerault/mineCETSAapp")*

*> library(mineCETSAapp)*

If you want now to use the app just type :

*> runCETSAapp()*
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
