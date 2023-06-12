# What is IMPRINTS.CETSA.app ?
This is a package to complete and provide a user friendly interface for the [IMRPRINTS.CETSA package](https://github.com/nkdailingyun/IMPRINTS.CETSA), more specifically for the analysis of the IMPRINTS-CETSA datasets.
With its Shiny app, the user can carry out the complete IMPRINTS-CETSA analysis at the protein and peptide level; visualize the bar plot or heatmap of proteins-of-interest; import new datasets or remove old ones at your end; get STRING network and perform GO analysis on your data; plot interactive networks of bar plots; start a PubMed search based on keywords such as protein names and save the results in Word files; etc.

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
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
