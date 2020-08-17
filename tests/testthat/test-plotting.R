context("plotting")

formatted_image <- SPIAT::formatted_image 

test_that("plot_average_expression() creates a plot", {
    
    p <- plot_average_expression(formatted_image, reference_marker="AMACR", target_marker="CD4", radii=c(30, 35, 40, 45, 50, 75, 100))
    
    expect_is(p, "ggplot")
})


test_that("plot_cell_categories() creates a plot", {
    
    phenotypes_of_interest <- c("AMACR", "CD3,CD4", "CD3,CD8")
    colour_vector <- c("darkgrey", "red", "blue")
    
    p <- plot_cell_categories(formatted_image, phenotypes_of_interest, colour_vector)
    
    expect_is(p, "ggplot")
})


test_that("plot_marker_level_heatmap() creates a plot", {
    
    p <- plot_marker_level_heatmap(formatted_image, num_splits = 100, "CD3")
    
    expect_is(p, "ggplot")
    
})


# test_that("plot_cell_marker_levels() creates a plot", {
# 
#     p <- plot_cell_marker_levels(formatted_image, return_data = FALSE)
# 
#     expect_is(p, "ggplot")
# 
# })


test_that("plot_distance_heatmap() creates a plot", {
    
    summary_distances <- calculate_summary_distances_between_phenotypes(formatted_image)
    
    p <- plot_distance_heatmap(summary_distances)
    
    expect_is(p, "ggplot")
    
})


# test_that("plot_cell_distances_violin() creates a plot", {
# 
#     distances <- calculate_all_distances_between_phenotypes(formatted_image,
#                                                             remove_other = TRUE,
#                                                             cell_phenotypes_of_interest = c("AMACR", "PDL1"))
# 
#     p <- plot_cell_distances_violin(distances)
# 
#     expect_is(p, "ggplot")
# 
# })


test_that("plot_composition_heatmap() creates a plot", {
    
    clusters <- identify_cell_clusters(formatted_image, phenotypes_of_interest = c("CD3,CD4", "CD3,CD8"),
                                       radius = 30)
    clusters_2 <- composition_of_clusters_and_communities(clusters, "Cluster")
    clusters_2 <- clusters_2[clusters_2$Total_number_of_cells >=5,]
    
    p <- plot_composition_heatmap(clusters_2, column_to_consider="Cluster")
    
    expect_is(p, "pheatmap")
    
})

