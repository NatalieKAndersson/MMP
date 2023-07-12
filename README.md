# MMP <img src="https://github.com/NatalieKAndersson/MMP/blob/main/Images/MMP_logo.PNG" align = "right" width="180"/>

The modified maximum parsimony (MMP) is an algorithm for phylogenetic reconstruction from multiregional sampling data that can incorporate information from SNP-array, WES, WGS, TDS etc. in unison or separately. It uses an event matrix along with information about each alteration's size in each sample as input and creates a phylogeny that makes sure the pigeonhole principle is not violated in any of them.

<a href="https://zenodo.org/badge/latestdoi/297145258"><img src="https://zenodo.org/badge/297145258.svg" alt="DOI"></a>

When using the algorithm, please cite: Rastegar, B., Andersson, N....

## Setting up MMP

MMP can be installed by running the following lines of code

```
library(devtools)
devtools::install_github('NatalieKAndersson/MMP')
library("MMP")
```


Before starting the analysis, make sure you have installed the dependencies of the algorithms.

```
library("readxl") #Needed to load the data from the xlsx file.
```
If they are not installed you can install them by using the following command.

```R
install.packages(c("readxl","xlsx","stringr","ape","phangorn","ggplot2",
                   "ggtree","ggimage","dplyr","RColorBrewer","ggridges","cowplot","dbscan"))
```


An alternative to installing the MMP package is to simply download the entire R script denoted "MMP.1.0.R" and double click on the script to open it in your R-environment. Then load all functions in the script by marking them and pressing “Run”.

If no error message has appeared, we are good to go!

## Usage
We are now ready to load some data and get going! At the top of the script you have to set the path i.e. the location of the files that will be analyzed.
```R
setwd("~/yourpath") #Should be where the DEVOLUTION-file is.
MMP_tree <- MMP(file="DEVOLUTION.xlsx",tumorname=x) #The output is the tree object.

- **File:** The filename of the output excel file from DEVOLUTION. This contains the event matrix, clustering and pie sizes. The algorithm uses the pie sizes etc. to minimize contradictions where “a larger pie comes after a smaller one” in the tree.
- **Tumorname:** The name of the tumor you are analyzing. In the function this is only used when the output file is saved in the end of the algorithm.
- The output of the MMP-function is the phylogenetic tree in the shape of a so called “phylo object” which can be plotted with ggplot.
```
