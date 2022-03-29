context("plotting")

test_that("plot_average_intensity() creates a plot", {
    
    p <-plot_average_intensity(simulated_image, reference_marker="Immune_marker3", 
                               target_marker="Immune_marker2", c(30, 35, 40, 45, 50, 75, 100))
    
    expect_is(p, "ggplot")
})

test_that("plot_cell_categories() creates a plot", {
    
    phenotypes_of_interest <- c("Tumour", "Immune2")
    colour_vector <- c("darkgrey", "blue")
    
    p <- plot_cell_categories(defined_image, phenotypes_of_interest, colour_vector,"Cell.Type")
    
    expect_is(p, "ggplot")
})

test_that("plot_cell_distances_violin() creates a plot", {
    
    pairwise_dist <- calculate_pairwise_distances_between_celltypes(
        defined_image, 
        cell_types_of_interest = c("Tumour","Immune1"), feature_colname = "Cell.Type")
    p <- plot_cell_distances_violin(pairwise_dist)
    
    expect_is(p, "ggplot")
})

test_that("plot_cell_percentages() creates a plot", {
    
    p_cells <- calculate_cell_proportions(simulated_image) 
    p <- plot_cell_percentages(p_cells)
    
    expect_is(p, "ggplot")
})

test_that("plot_cell_marker_levels() creates a plot", {
    
    p <- plot_cell_marker_levels(simulated_image, "Immune_marker1")
    
    expect_is(p, "ggplot")
    
})

test_that("plot_marker_level_heatmap() creates a plot", {
    
    p <- plot_marker_level_heatmap(simulated_image, num_splits = 100, "Tumour_marker")
    
    expect_is(p, "ggplot")
    
})

test_that("plot_distance_heatmap() creates a plot", {
    pairwise_dist <- calculate_pairwise_distances_between_celltypes(defined_image, 
        cell_types_of_interest = c("Tumour","Immune1"), feature_colname = "Cell.Type")
    
    summary_distances <- calculate_summary_distances_between_celltypes(pairwise_dist)
    p <- plot_distance_heatmap(summary_distances)
    
    expect_is(p, "ggplot")
    
})

test_that("marker_intensity_boxplot() creates a plot", {
    
    p <- marker_intensity_boxplot(simulated_image, "Immune_marker1")
    
    expect_is(p, "ggplot")
    
})

test_that("marker_prediction_plot() creates a plot", {
    predicted_result <- predict_phenotypes(sce_object = simulated_image, 
                                           thresholds = NULL,tumour_marker = "Tumour_marker",
                                           baseline_markers = c("Immune_marker1", "Immune_marker2","Immune_marker3", "Immune_marker4"), 
                                           reference_phenotypes = TRUE)
    p <- marker_prediction_plot(predicted_result, marker = "Tumour_marker")
    expect_is(p, "gtable")    
})

test_that("marker_surface_plot() creates a plot", {
    
    p <- marker_surface_plot(simulated_image, num_splits=15, marker="Immune_marker1")
    
    expect_is(p, "plotly")
    
})

test_that("marker_surface_plot_stack() creates a plot", {
    
    p <- marker_surface_plot_stack(simulated_image, num_splits=15, 
                                   markers=c("Tumour_marker", "Immune_marker4"))
    
    expect_is(p, "plotly")
})

test_that("dimensionality_reduction_plot() creates a plot", {
    
    p <- dimensionality_reduction_plot(simulated_image, plot_type = "TSNE", 
                                       feature_colname = "Phenotype")
    
    expect_is(p, "ggplot")
})

