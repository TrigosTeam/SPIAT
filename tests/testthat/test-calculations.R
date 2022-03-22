context("calculations")

test_that("average_marker_intensity_within_radius() works", {
    res <- 0.5995
    out <- average_marker_intensity_within_radius(simulated_image,
               reference_marker ="Immune_marker3", target_marker = "Immune_marker2",
               radius=30)
    expect_equal(out, res, tolerance=1e-4)
})

test_that("average_minimum_distance() works", {
    res <- 17.0134
    out <- average_minimum_distance(simulated_image)
    expect_equal(out, res, tolerance=1e-4)
})

test_that("average_nearest_neighbor_index() works", {
    res <- list()
    res$pattern <- "Clustered"
    res$`p-value` <- 7.11006e-32
    out <- average_nearest_neighbor_index(defined_image, "Tumour", "Cell.Type")
    expect_equivalent(out, res)
})

test_that("average_percentage_of_cells_within_radius() works", {
    res <- 17.4571
    out <- average_percentage_of_cells_within_radius(defined_image, "Tumour", 
                                                     "Immune3", radius = 100, "Cell.Type")
    expect_equal(out, res, tolerance=1e-4)
})

test_that("calculate_cell_proportions() works", {
    res <- data.frame(row.names = c(4L, 5L, 2L, 1L,3L),
                          Cell_type = factor(c("OTHER", "Tumour_marker", "Immune_marker1,Immune_marker2,Immune_marker4", 
                                               "Immune_marker1,Immune_marker2", "Immune_marker1,Immune_marker3")), 
                          Number_of_celltype = c(2986, 819, 630, 338, 178),
                          Proportion = c(0.60311048273076151016,  0.16542112704504141618, 
                                         0.12724702080387800818, 0.06826903655827105954, 0.03595233286204806838), 
                          Percentage = c(60.31104827307615323662, 16.54211270450414161814, 
                                         12.72470208038780015158, 6.82690365582710612102, 3.59523328620480686624),
                          Proportion_name = rep("/Total", 5),
                          stringsAsFactors = FALSE)
    
    p_cells <- calculate_cell_proportions(sce_object = simulated_image)
    
    expect_equal(p_cells, res)
})

test_that("calculate_distances_between_cell_types() works", {
    res <- data.frame(row.names = c(2L, 3L, 4L),
                      Cell1 = c("Cell_25", "Cell_30", "Cell_48"), 
                      Cell2 = c("Cell_15", "Cell_15", "Cell_15"),
                      Distance = c(119.3982738563,  61.7785923437, 279.8494352805), 
                      Pair = rep("Immune1_Immune1", 3))
    
    dists <- calculate_distances_between_cell_types(defined_image, 
                                                    cell_types_of_interest = c("Tumour","Immune1"), 
                                                    feature_colname = "Cell.Type")
    out <- dists[1:3, ]
    out$Cell1 <- as.character(out$Cell1)
    out$Cell2 <- as.character(out$Cell2)
    expect_equal(out, res, tolerance = 1e-6)
})

test_that("calculate_summary_distances_between_cell_types() works", {
    
    res <- data.frame(row.names = c(1L, 2L),
                      Pair = c("Immune1_Immune1", "Immune1_Tumour"),
                      Mean = c(1164.7096, 1013.3697), 
                      Min = c(10.84056, 13.59204),
                      Max = c(2729.120, 2708.343),
                      Median = c(1191.3645, 1004.6579),
                      Std.Dev = c(552.0154, 413.7815),
                      Reference = c("Immune1", "Immune1"),
                      Target = c("Immune1", "Tumour"),
                      stringsAsFactors = FALSE)
    
    summary_distances <- calculate_summary_distances_between_cell_types(SPIAT::defined_image,
                            feature_colname = "Cell.Type", all_combinations = FALSE,
                            cell_types_of_interest = c("Tumour","Immune1"))
    out <- summary_distances[c(1,2),]
    
    expect_equal(out, res, tolerance=1e-4)
})

test_that("calculate_entropy() works", {
    out <- calculate_entropy(defined_image, cell_types_of_interest = c("Immune1","Immune2"), 
                             feature_colname = "Cell.Type")
    expect_equal(0.9294873, out, tolerance=1e-4)
})

test_that("calculate_minimum_distances_between_cell_types() works", {
    res <- data.frame(row.names = c(2L, 3L),
                      Reference = c("Immune1", "Tumour"),
                      Nearest = c("Tumour", "Immune1"),
                      Mean = c(85.84773, 187.52474), 
                      Std.Dev = c(40.72454, 75.13630),
                      Median = c(80.80592, 191.09109),
                      stringsAsFactors = FALSE)
    out <- calculate_minimum_distances_between_cell_types(defined_image,
                                                         feature_colname = "Cell.Type", 
                                                         cell_types_of_interest = c("Tumour","Immune1"))
    expect_equal(res, out, tolerance=1e-4)
})

test_that("calculate_minimum_distances() works", {
    res <- data.frame(row.names = c(2L, 3L, 4L),
                      RefCell = c("Cell_15", "Cell_25", "Cell_30"),
                      RefType = c("Immune1", "Immune1", "Immune1"),
                      NearestCell = c("Cell_32", "Cell_27", "Cell_32"),
                      NearestType = c("Tumour", "Tumour", "Tumour"),
                      Dist = c(17.18740, 44.79503, 78.52918))
    min_dists <- calculate_minimum_distances(defined_image, 
                                             cell_types_of_interest = c("Tumour","Immune1"),
                                             feature_colname = "Cell.Type")
    out <- min_dists[1:3,]
    expect_equal(res, out, tolerance=1e-4)
})

test_that("marker_permutation() works", {
    
    res <- data.frame(row.names = c("Tumour_marker", "Immune_marker1", "Immune_marker2"),
                      Observed_cell_number = c(819, 0, 0),
                      Percentage_of_iterations_where_present = c(100, 100, 100),
                      Average_bootstrap_cell_number = c(423.41, 643.25, 517.77),
                      Enrichment.p = c(0.01, 1.00, 1.00), 
                      Depletion.p = c(1.00, 0.01, 0.01))
    
    sig <- marker_permutation(simulated_image, num_iter = 100)
    out <- sig[1:3,]
    
    expect_equal(out, res, tolerance=1e-2)
    
})

test_that("mixing_score_summary() works", {
    
    res <- data.frame(row.names = c(2L),
                      Reference = "Tumour",
                      Target = "Immune1",
                      Number_of_reference_cells = 819,
                      Number_of_target_cells = 338,
                      Reference_target_interaction = 80,
                      Reference_reference_interaction = 5026,
                      Mixing_score = 0.0159,
                      Normalised_mixing_score = 0.077)
    
    out <- mixing_score_summary(defined_image, reference_celltype = "Tumour", target_celltype="Immune1",
                                radius = 50, feature_colname = "Cell.Type")
    
    expect_equal(out, res, tolerance=1e-2)
    
})

test_that("number_of_cells_within_radius() works", {
    
    res <- data.frame(row.names = c("Cell_19", "Cell_27"),
                      Cell.X.Position = c(11.82661, 144.71614),
                      Cell.Y.Position = c(145.4592, 172.2780),
                      Immune1 = c(0,1))
    
    n_in_radius <- number_of_cells_within_radius(defined_image, 
                                                reference_celltype = "Tumour", 
                                                target_celltype="Immune1", radius = 50, 
                                                feature_colname = "Cell.Type")
    out <- n_in_radius$Tumour[1:2,]
    
    expect_is(n_in_radius, "list")
    expect_equal(out, res, tolerance=1e-2)
})

test_that("compute_gradient() works", {
    
    gradient_positions <- c(30, 50, 100)
    gradient_entropy <- compute_gradient(defined_image, radii = gradient_positions, 
                                         FUN = calculate_entropy,  
                                         cell_types_of_interest = c("Immune1","Immune2"),
                                         feature_colname = "Cell.Type")

    out <- gradient_entropy[[1]]
    
    expect_is(gradient_entropy, "list")
    expect_length(gradient_entropy, 3)
    expect_equal(dim(out), c(338,13))
})

test_that("compute_gradient() works", {
    
    res <- data.frame(row.names = c(1L, 2L),
                      Celltype1 = c("Tumour", "Immune3"),
                      Celltype2 = c("Immune3", "Tumour"),
                      Pos_50 = c(0.5974227,0.8053223),
                      Pos_100 = c(0.7667379, .9667266))
    
    gradient_pos <- seq(50, 500, 50)
    gradient_results <- entropy_gradient_aggregated(defined_image, 
                                                    cell_types_of_interest = c("Tumour","Immune3"),
                                                    feature_colname = "Cell.Type", radii = gradient_pos)
    
    out <- gradient_results$gradient_df[, 1:4]
    
    expect_is(gradient_results, "list")
    expect_equal(gradient_results$peak, 10)
    expect_equal(res, out, tolerance = 1e-4)
})

test_that("measure_association_to_cell_properties() works", {
    t <- measure_association_to_cell_properties(image_no_markers,
                                                celltypes = c("Tumour", "Immune2"),
                                                feature_colname="Cell.Type",
                                                property = "Cell.Size",
                                                method = "t")
    expect_is(t, "htest")
    expect_equal(unname(t$statistic), 27.4, tolerance = 0.1)
    
    p <- measure_association_to_cell_properties(image_no_markers,
                                                celltypes = c("Tumour", "Immune2"),
                                                feature_colname="Cell.Type",
                                                property = "Cell.Size",
                                                method = "box")
    expect_is(p, "ggplot")
})
