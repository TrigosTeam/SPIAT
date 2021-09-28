context("calculations")

test_that("calculate_cell_proportions() works", {
    
    res <- data.frame(row.names = c(1L, 5L, 3L, 4L, 2L, 6L),
                          Cell_type = factor(c("AMACR", "OTHER", "CD3,CD4", "CD3,CD8", "CD3", "PDL-1")), 
                          Number_of_celltype = c(4446, 3299, 513, 138, 19, 4),
                          Proportion = c(0.52809122223542, 
                                     0.391851763867443, 0.0609336025656254, 0.0163914954270103, 0.00225680009502316, 
                                     0.00047511580947856), 
                          Percentage = c(52.809122223542, 39.1851763867443, 
                                     6.09336025656254, 1.63914954270103, 0.225680009502316, 0.047511580947856),
                          Proportion_name = c("AMACR/Total", "OTHER/Total", "CD3,CD4/Total", "CD3,CD8/Total", "CD3/Total", "PDL-1/Total"),
                          stringsAsFactors = FALSE)
    
    p_cells <- calculate_cell_proportions(sce_object = formatted_image)
    
    expect_equal(p_cells, res)
    
})


test_that("compute_mixing_score() works", {
    
    score <- compute_mixing_score(formatted_image, reference_marker = "CD4", target_marker = "PDL-1")
    
    expect_equal(score, 0.009708738)
    
})


test_that("average_marker_intensity_within_radius() works", {
    
    ave_marker_intensity <- average_marker_intensity_within_radius(formatted_image,
                                                   reference_marker ="AMACR",
                                                   target_marker = "CD3",
                                                   radius=30)
    
    expect_equal(ave_marker_intensity, 3.717, tolerance=1e-3)
    
})


test_that("average_minimum_distance() works", {
    
    ave_min_dist <- average_minimum_distance(formatted_image)
    
    expect_equal(ave_min_dist, 15.924, tolerance=1e-3)
    
})


test_that("calculate_summary_distances_between_cell_types() works", {
    
    res <- data.frame(Reference = c("OTHER", "OTHER", "OTHER"),
                      Nearest = c("AMACR", "CD3,CD4", "CD3,CD8"),
                      Mean = c(130.418076521445, 70.5496420825272, 134.633340256936), 
                      Std.Dev = c(123.531481047835, 61.6136499595389, 124.074239935635),
                      Median = c(85.8661749468322, 52.3545604508337, 97.0463806640928))
    
    summary_distances <- calculate_summary_distances_between_cell_types(formatted_image)
    rownames(summary_distances) <- NULL
    
    expect_equal(summary_distances[1:3, ], res, tolerance=1e-3)
    
})


test_that("calculate_all_distances_between_cell_types() works", {
    
    cells <- c("Cell_174", "Cell_264", "Cell_305")
    res <- data.frame(Cell1 = factor(cells, levels=cells),  
                      Cell2 = factor(c(rep("Cell_78", 3))),
                      Distance = c(1392.07040, 40.26164, 993.83399), 
                      Pair = rep("CD3,CD4_CD3,CD4", 3))
    
    distances <- calculate_all_distances_between_cell_types(formatted_image,
                                                            cell_types_of_interest = c("CD3,CD4", "CD3,CD8"),
                                                            column="Phenotype")
    distances <- distances[1:3, ]
    distances <- droplevels(distances)
    rownames(distances) <- NULL
    
    expect_equal(distances, res, tolerance=1e-3)
    
})


test_that("identify_cell_clusters() works", {
    
    cells <- c("Cell_78", "Cell_174", "Cell_195", "Cell_237")
    res <- data.frame(row.names = c("Cell_195", "Cell_237", "Cell_320", "Cell_441"),
                Cell.ID = c("Cell_195", "Cell_237", "Cell_320", "Cell_441"), 
                Phenotype = c("CD3,CD8", "CD3,CD8", "CD3,CD8", "CD3,CD8"), 
                Cell.X.Position = c(1145, 1119, 550, 856), 
                Cell.Y.Position = c(33, 40, 73, 93))
    
    clusters <- identify_cell_clusters(formatted_image, cell_types_of_interest = c("CD3,CD8"),
                                       radius = 30, column = "Phenotype")
    
    expect_equal(clusters[1:4, 1:4], res)
    
})

test_that("identify_cell_communities() works", {
    
    res <- data.frame(row.names = 2:5,
                      Cell.ID = c("Cell_2", "Cell_3", "Cell_4", "Cell_5"), 
                      Phenotype = c("AMACR", "AMACR", "AMACR", "AMACR"), 
                      Cell.X.Position = c(171, 184, 201, 219), 
                      Cell.Y.Position = c(22, 38, 52, 63))
    
    communities <- identify_cell_communities(formatted_image, radius=100)
    
    expect_equal(communities[1:4, 1:4], res)
    
})


test_that("marker_permutation() works", {
    
    res <- data.frame(row.names = c("CD3", "PDL-1", "CD4", "CD8"),
                      Observed_cell_number = c(19, 4, 0, 0),
                      Percentage_of_iterations_where_present = c(100, 88, 100, 100),
                      Average_bootstrap_cell_number = c(292.52, 1.71, 216.97, 55.81),
                      Enrichment.p = c(1, 0.05, 1, 1), 
                      Depletion.p = c(0.01, 1, 0.01, 0.01))
    
    sig <- marker_permutation(formatted_image, num_iter = 100)
    sig <- sig[1:4,]
    
    expect_equal(sig, res, tolerance=0.9)
    
})


test_that("average_percentage_of_cells_within_radius() works", {
    
    percent_cells <- average_percentage_of_cells_within_radius(formatted_image, reference_phenotypes = "PDL-1", target_phenotypes = "AMACR", radius=100, column="Phenotype")
    
    expect_equal(percent_cells, 10.52296, tolerance=1e-7)
    
})
