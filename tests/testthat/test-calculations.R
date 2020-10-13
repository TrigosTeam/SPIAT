context("calculations")

test_that("calculate_cell_proportions() works", {
    
    res <- data.frame(row.names = c(1L, 5L, 3L, 4L, 2L, 6L),
                          Cell_type = factor(c("AMACR", "OTHER", "CD3,CD4", "CD3,CD8", "CD3", "PDL-1")), 
                          Number_of_celltype = c(4446, 3299, 513, 138, 19, 4),
                          Reference = rep("Total", 6),
                          Number_of_reference = rep(8419, 6),
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


test_that("calculate_summary_distances_between_phenotypes() works", {
    
    res <- data.frame(Target = c("OTHER", "OTHER", "OTHER"),
                      Nearest = c("AMACR", "CD3,CD4", "CD3,CD8"),
                      Mean = c(130.418076521445, 70.5496420825272, 134.633340256936), 
                      Std.Dev = c(123.531481047835, 61.6136499595389, 124.074239935635),
                      Median = c(85.8661749468322, 52.3545604508337, 97.0463806640928))
    
    summary_distances <- calculate_summary_distances_between_phenotypes(formatted_image)
    rownames(summary_distances) <- NULL
    
    expect_equal(summary_distances[1:3, ], res, tolerance=1e-3)
    
})


test_that("calculate_all_distances_between_phenotypes() works", {
    
    cells <- c("Cell_174", "Cell_264", "Cell_305")
    res <- data.frame(Var1 = factor(cells, levels=cells),  
                      Var2 = factor(c(rep("Cell_78", 3))),
                      Distance = c(1392.07040, 40.26164, 993.83399), 
                      Pair = rep("CD3,CD4_CD3,CD4", 3))
    
    distances <- calculate_all_distances_between_phenotypes(formatted_image,
                                                            remove_other = TRUE,
                                                            cell_phenotypes_of_interest = c("CD3,CD4", "CD3,CD8"))
    distances <- distances[1:3, ]
    distances <- droplevels(distances)
    rownames(distances) <- NULL
    
    expect_equal(distances, res, tolerance=1e-3)
    
})


test_that("identify_cell_clusters() works", {
    
    cells <- c("Cell_78", "Cell_174", "Cell_195", "Cell_237")
    res <- data.frame(row.names = cells,
                      Cell.ID = cells,
                      Phenotype = c("CD3,CD4", "CD3,CD4", "CD3,CD8", "CD3,CD8"),
                      Cell.X.Position = c(1079, 2471, 1145, 1119),
                      Cell.Y.Position = c(15, 29, 33, 40),
                      Cell.Area = c(129, 448, 126, 223),
                      Nucleus.Area = c(87, 144, 43, 22),
                      Nucleus.Compactness = c(0.84, 0.42, 0.8, 1.08),
                      Nucleus.Axis.Ratio = c(1, 2.16, 1.14, 1),
                      Cell.Axis.Ratio = c(1.17, 1.41, 1.83, 2.03),
                      DAPI = c(14.7, 12.8, 12.7, 11),
                      CD3 = c(3.67, 2.39, 1.3, 1.38),
                      `PDL-1` = c(0.357, 0.2, 0.114, 0.23),
                      CD4 = c(3.74, 1.22, 1.9, 1.3),
                      CD8 = c(0.319, 0.05, 9.982, 4.254),
                      AMACR = c(0.034, 0.069, 0.514, 0.061),
                      Cluster = c("Cluster_1", "Cluster_NA", "Cluster_1", "Cluster_1"), check.names = FALSE)
    
    clusters <- identify_cell_clusters(formatted_image, phenotypes_of_interest = c("CD3,CD4", "CD3,CD8"),
                                       radius = 30)
    
    expect_equal(clusters[1:4,], res, tolerance = 0.002)
    
})


test_that("identify_cell_communities() works", {
    
    res <- data.frame(row.names = 2:5,
                      Cell.ID = c("Cell_2", "Cell_3", "Cell_4", "Cell_5"), 
                      Phenotype = c("AMACR", "AMACR", "AMACR", "AMACR"), 
                      Cell.X.Position = c(171L, 184L, 201L, 219L), 
                      Cell.Y.Position = c(22L, 38L, 52L, 63L), 
                      Cell.Area = c(464L, 553L, 462L, 876L), 
                      Nucleus.Area = c(177L, 212L, 239L, 451L), 
                      Nucleus.Compactness = c(0.54, 0.51, 0.53, 0.53), 
                      Nucleus.Axis.Ratio = c(1.84, 1.92, 1.47, 1.19), 
                      Cell.Axis.Ratio = c(1, 1.21, 1.09, 1.34), 
                      DAPI = c(17.588, 21.262, 18.951, 18.631), 
                      CD3 = c(0.229, 0.206, 0.226, 0.212), 
                      `PDL-1` = c(0.117, 0.103, 0.126, 0.091), 
                      CD4 = c(1.188, 0.924, 1.367, 1.266), 
                      CD8 = c(0.087, 0.053, 0.037, 0.046), 
                      AMACR = c(2.074, 1.521, 2.462, 1.968), 
                      Community = c("Community_1", "Community_1", "Community_1", "Community_1"), check.names = FALSE)
    
    communities <- identify_cell_communities(formatted_image, radius=100)
    
    expect_equal(communities[1:4,], res, tolerance=0.002)
    
})


test_that("marker_permutation() works", {
    
    res <- data.frame(row.names = c("CD3", "PDL-1", "CD4", "CD8"),
                      Percentage_of_occurrence = c(100, 90, 100, 100),
                      Observed_cell_number = c(19, 4, 0, 0),
                      Average_bootstrap_cell_number = c(291.78, 1.75, 216.96, 54.72),
                      Enrichment.p = c(1, 0.06, 1, 1), 
                      Depletion.p = c(0.01, 1, 0.01, 0.01))
    
    sig <- marker_permutation(formatted_image, num_iter = 100)
    sig <- sig[1:4,]
    
    expect_equal(sig, res, tolerance=0.9)
    
})


test_that("percentage_of_cells_within_radius() works", {
    
    res <- setNames(c(5.357143, 36.734694, 0.00000, 0.00000 ), c("Cell_664","Cell_2760", "Cell_2992", "Cell_3147"))
    
    percent_cells <- percentage_of_cells_within_radius(formatted_image, reference_phenotypes = "PDL-1", target_phenotypes = "AMACR", radius=100)
    
    expect_equal(percent_cells, res, tolerance=1e-7)
    
})
