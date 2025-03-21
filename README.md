# SPIAT - Spatial Image Analysis of Tissues

SPIAT (Spatial Image Analysis of Tissues) is an R package with a suite of data processing, quality control, visualization, data handling and data analysis tools. SPIAT is compatible with data generated from single-cell spatial proteomics platforms (e.g. OPAL, CODEX, MIBI). SPIAT reads spatial data in the form of X and Y coordinates of cells, marker intensities and cell phenotypes.

SPIAT includes six analysis modules that allow visualization, calculation of cell colocalization, categorization of the immune microenvironment relative to tumor areas, analysis of cellular neighborhoods, and the quantification of spatial heterogeneity, providing a comprehensive toolkit for spatial data analysis.


## Installation

To install this package, start R and enter:
```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("SPIAT")
```
You can also install the latest development version from Github.

```r
if (!requireNamespace("devtools", quietly = TRUE)) {
      install.packages("devtools")}
devtools::install_github("TrigosTeam/SPIAT")
```
The estimated installation time on a Windows 10 (64-bit) system is 1.5 minutes.

## Vignette 

The vignette with an overview of the package can be accessed from the top Menu under Articles or by clicking [here](https://trigosteam.github.io/SPIAT/articles/SPIAT.html).

## Issues

Please open an issue on our Github page (https://github.com/TrigosTeam/SPIAT/issues) if you have any questions.

## Authors of the package
Yuzhou Feng, Tianpei Yang, Volkan Ozcoban, Mabel Li and John Zhu developed the package, including developing algorithms, writing code and designing the package. Yuzhou Feng and Maria Doyle did the package cleaning and wrote the tutorial. Anna Trigos conceived, supervised the work, developed algorithms and designed the package. 

## Paper reference
Please check our latest paper for more information!

Yuzhou Feng et al, Spatial analysis with SPIAT and spaSim to characterize and simulate tissue microenvironments, Nature Communications (2023). DOI: 10.1038/s41467-023-37822-0
