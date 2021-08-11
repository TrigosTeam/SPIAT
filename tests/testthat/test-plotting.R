context("plotting")

test_that("plot_average_intensity() creates a plot", {
    
    p <- plot_average_intensity(formatted_image, reference_marker="AMACR", target_marker="CD4", radii=c(30, 35, 40, 45, 50, 75, 100))
    
    expect_is(p, "ggplot")
})


test_that("plot_cell_categories() creates a plot", {
    
    phenotypes_of_interest <- c("AMACR", "CD3,CD4", "CD3,CD8")
    colour_vector <- c("darkgrey", "red", "blue")
    
    p <- plot_cell_categories(formatted_image, phenotypes_of_interest, colour_vector)
    
    expect_is(p, "ggplot")
})

test_that("plot_cell_percentages() creates a plot", {
    
    p_cells <- calculate_cell_proportions(sce_object = formatted_image)
    
    p <- plot_cell_percentages(p_cells)
    
    expect_is(p, "ggplot")
})


test_that("plot_cell_marker_levels() creates a plot", {
    
    p <- plot_cell_marker_levels(formatted_image, "CD3")
    
    expect_is(p, "ggplot")
    
})


test_that("plot_marker_level_heatmap() creates a plot", {
    
    p <- plot_marker_level_heatmap(formatted_image, num_splits = 100, "CD3")
    
    expect_is(p, "ggplot")
    
})


test_that("plot_distance_heatmap() creates a plot", {
    
    summary_distances <- calculate_summary_distances_between_cell_types(formatted_image)
    
    p <- plot_distance_heatmap(summary_distances)
    
    expect_is(p, "ggplot")
    
})


test_that("plot_composition_heatmap() creates a plot", {
    
    clusters <- identify_cell_clusters(formatted_image, cell_types_of_interest = c("CD3,CD4", "CD3,CD8"),
                                       radius = 30, column = "Phenotype")
    clusters_2 <- composition_of_clusters_and_communities(clusters, type_of_aggregate = "Cluster", column = "Phenotype")
    clusters_2 <- clusters_2[clusters_2$Total_number_of_cells >=5,]
    
    p <- plot_composition_heatmap(clusters_2, type_of_aggregate = "Cluster", column="Phenotype")
    
    expect_is(p, "Heatmap")
    
})

# Comment out until update with Jojo's new function
# test_that("identify_bordering_cells() creates a plot", {
#     
#     p <- identify_bordering_cells(formatted_image, reference_marker = "AMACR",
#                                   rm_noise_radius = 50, radius = 100, lower_bound = 0.05,
#                                   upper_bound=0.7)
#     
#     expect_is(p, "ggplot")
#     
# })


test_that("marker_intensity_boxplot() creates a plot", {
    
    p <- marker_intensity_boxplot(formatted_image, "CD3")
    
    expect_is(p, "ggplot")
    
})


#test_that("marker_prediction_plot() creates a plot", {
#    
#    predicted_image <- predict_phenotypes(formatted_image,
#                                          thresholds = NULL,
#                                          tumour_marker = "AMACR",
#                                          baseline_markers = c("CD3", "CD4", "CD8"),
#                                          reference_phenotypes = FALSE)
#    
#    p <- marker_prediction_plot(predicted_image, marker="AMACR")
#    
#    expect_is(p, "gtable")
#    
#})


test_that("marker_surface_plot() creates a plot", {
    
    p <- marker_surface_plot(formatted_image, num_splits=15, marker="CD3")
    
    expect_is(p, "plotly")
    
})


test_that("marker_surface_plot_stack() creates a plot", {
    
    p <- marker_surface_plot_stack(formatted_image, num_splits=10, markers=c("CD4", "AMACR"))
    
    expect_is(p, "plotly")
    
})


test_that("measure_association_to_cell_properties() creates a plot", {
    
    p <- measure_association_to_cell_properties(formatted_image, phenotypes = c("CD3,CD4", "CD3,CD8"))
    
    expect_is(p, "ggplot")
    
})
