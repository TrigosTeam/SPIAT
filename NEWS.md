# SPIAT 0.99.1

SIGNIFICANT USER-VISIBLE CHANGES

* Updated the main object class from `SingleCellExperiment` to `SpatialExperiment`.
* Replace "sce" with "spe" in function names.
* Deleted "Visium" format option in `format_image_to_spe()`.
* `image_splitter()` now returns a list of spe objects.

# SPIAT 0.99.0

Version submitted to Nature Communications. Access the release v0.99.0 [here](https://github.com/TrigosTeam/SPIAT/releases/tag/v0.99.0).

SIGNIFICANT USER-VISIBLE CHANGES

* Systematically renamed parameter names. For example, the parameter for selecting a column of interest is changed from "column" to "feature_colname".
* Added a function `entropy_gradient_aggregated()` for calculating cell colocalisation. 
* Added a function `dimensionality_reduction_plot()` for visualisation and quality control.
* Added examples to documentation. The example datasets were simulated by [**spaSim**](https://github.com/TrigosTeam/spaSim/releases/tag/v0.99.1).

# SPIAT 0.4

SIGNIFICANT USER-VISIBLE CHANGES

* Added a new parameter "column" to most functions to select a column of interest.
* Users can define cell types based on certain marker combinations or cell properties.
* Improved tumour border detection (`identify_bordering_cells()`).
* Added functions to define tumour structure and characterise immune populations in different tissue regions.
* Added functions to characterise spatial heterogeneity.
* Added functions for cell colocalisation, including mixing score, normalised mixing score, cross K function.

# SPIAT 0.2

* Version for the [paper](https://www.biorxiv.org/content/10.1101/2020.05.28.122614v1) submitted to BioRxiv.
