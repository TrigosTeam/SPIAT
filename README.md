# SPIAT - Spatial Image Analysis of Tissues

SPIAT (Spatial Image Analysis of Tissues) is an R package with a suite of data processing, quality control, visualization, data handling and data analysis tools. SPIAT is compatible with data generated from single-cell spatial proteomics platforms (e.g. OPAL, CODEX, MIBI). SPIAT reads spatial data in the form of X and Y coordinates of cells, marker intensities and cell phenotypes.

SPIAT includes six analysis modules that allow visualization, calculation of cell colocalization, categorization of the immune microenvironment relative to tumor areas, analysis of cellular neighborhoods, and the quantification of spatial heterogeneity, providing a comprehensive toolkit for spatial data analysis.


## Installation

```r
# install.packages('devtools')
devtools::install_github("TrigosTeam/SPIAT")
```

## Vignette 

The vignette with an overview of the package can be accessed by clicking [here](https://trigosteam.github.io/SPIAT/articles/introduction.html)

## Authors of the package
Yuzhou Feng, Tianpei Yang, Volkan Ozcoban, Mabel Li and John Zhu developed the package, including developing algorithms, writing code and desgining the package. Anna Trigos conceived, supervised the work, developed algorithms and desgined the package.

## Paper reference
Our latest paper will be submitted shortly! However, you can have a look at our previous SPIAT pre-print here:
https://www.biorxiv.org/content/10.1101/2020.05.28.122614v1
