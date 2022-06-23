# SPIAT - Spatial Image Analysis of Tissues

SPIAT (Spatial Image Analysis of Tissues) is an R package with a suite of data processing, quality control, visualization, data handling and data analysis tools. SPIAT is compatible with data generated from single-cell spatial proteomics platforms (e.g. OPAL, CODEX, MIBI). SPIAT reads spatial data in the form of X and Y coordinates of cells, marker intensities and cell phenotypes.

SPIAT includes six analysis modules that allow visualization, calculation of cell colocalization, categorization of the immune microenvironment relative to tumor areas, analysis of cellular neighborhoods, and the quantification of spatial heterogeneity, providing a comprehensive toolkit for spatial data analysis.


## Installation

To install this package, start R and enter:
```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("SPIAT")
```

The estimated installation time on a Windows 10 (64-bit) system is 1.5 minutes.

## Vignette 

The vignette with an overview of the package can be accessed from the top Menu under Articles or by clicking [here](https://trigosteam.github.io/SPIAT/articles/introduction.html).

## Authors of the package
Yuzhou Feng, Tianpei Yang, Volkan Ozcoban, Mabel Li and John Zhu developed the package, including developing algorithms, writing code and designing the package. Yuzhou Feng and Maria Doyle did the package cleaning and wrote the tutorial. Anna Trigos conceived, supervised the work, developed algorithms and desgined the package. 

## Paper reference
Our latest paper is currently under review! However, you can have a look at our previous SPIAT pre-print here:
https://www.biorxiv.org/content/10.1101/2020.05.28.122614v1

Note that the version used in the manuscript under review is v0.99.0. The release can be accessed [here](https://github.com/TrigosTeam/SPIAT/releases/tag/v0.99.0).

