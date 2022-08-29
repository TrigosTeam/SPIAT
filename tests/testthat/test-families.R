context("families")

test_that("functions in cross K family work", {
    ## calculate_cross_functions()
    res <- data.frame(row.names = c(10L, 11L, 12L),
                      r = c(9.09090909090, 10.10101010101, 11.11111), 
                      theo = c(259.6357564950, 320.53797098151,387.850944887),
                      border = c(0.0000000000,7.75578194807,77.5578194807),
                      stringsAsFactors = FALSE)
    df_cross <- calculate_cross_functions(
        defined_image, method = "Kcross",
        cell_types_of_interest = c("Tumour","Immune3"),
        feature_colname ="Cell.Type", dist = 100, plot_results = FALSE)
     # test if the data strcure is "fv"
    expect_s3_class(df_cross, "fv")   
    
    out <- data.frame(df_cross[10:12,])
    expect_equal(res, out, tolerance = 1e-4) # test if the results are the same
    
    ## AUC_of_cross_function()
    res <- 0.07372796865617478601
    out <- AUC_of_cross_function(df_cross)
    expect_equal(res, out)
    
    ## crossing_of_crossK(df_cross)
    res <- 0.5
    out <- crossing_of_crossK(df_cross)
    expect_equal(res, out)
})

test_that("functions in tumour structure family work", {
    ## identify_bordering_cells()
    res <- c("Inside", "Inside", "Border", "Inside")
    spe_border <- identify_bordering_cells(
        defined_image, reference_cell = "Tumour", feature_colname = "Cell.Type",
        n_to_exclude = 10)
    # test if the data strcure is "SpatialExperiment"
    expect_s4_class(spe_border, "SpatialExperiment") 
    
    out <- spe_border$Region[3776:3779]
    expect_equal(res, out) # test if the results are the same
    
    ## R_BT()
    # test if R_BT returns a number
    n <- R_BT(defined_image, cell_type_of_interest = "Tumour", "Cell.Type")
    expect_is(n, "numeric")
    
    ## calculate_distance_to_tumour_margin()
    res <- c(83.54843614991348488275, 63.04062766766428183018, 
             0.00000000000000000000, 29.59984593965495847101)
    res1 <- c("Non-border", "Non-border", "Border", "Non-border")
    spe_dist <- calculate_distance_to_tumour_margin(spe_border)
    # test if the data strcure is "SpatialExperiment"
    expect_s4_class(spe_dist, "SpatialExperiment") 
    out <- spe_dist$Distance.To.Border[3776:3779]
    out1 <- spe_dist$region2[3776:3779]
    # test if the results are the same
    expect_equal(res, out) 
    expect_equal(res1, out1)
    
    ## define_structure()
    res <- c("Internal.margin", "Internal.margin", "Border", "Internal.margin")
    spe_structure <- define_structure(
        spe_dist, names_of_immune_cells = c("Immune1","Immune2","Immune3"),
        feature_colname = "Cell.Type", n_margin_layers = 5)
    
    # test if the data strcure is "SpatialExperiment"
    expect_s4_class(spe_structure, "SpatialExperiment") 
    
    out <- spe_structure$Structure[3776:3779]
    expect_equal(res, out)  # test if the results are the same
    
    ## summarise the cell proportions in each tumour structure
    res <- data.frame(
        Cell.Type = c("Immune1", "Immune3", "Immune1", "Immune3", 
                      "Immune1", "Immune3", "All_cells_of_interest"),
        Relative_to = c("All_cells_in_the_structure", 
                        "All_cells_in_the_structure", 
                        "All_cells_of_interest_in_the_structure", 
                        "All_cells_of_interest_in_the_structure", 
                        "The_same_cell_type_in_the_whole_image", 
                        "The_same_cell_type_in_the_whole_image", 
                        "All_cells_in_the_structure"),
        P.Infiltrated.Immune=as.numeric(c("0", "0.12576687", "0", "1", 
                                          "0", "0.06507937", "0.14385965")),
        P.Internal.Margin.Immune = as.numeric(c("0", "0.08071749", "0", "1", 
                                                "0", "0.05714286", "0.08780488")),
        P.External.Margin.Immune = as.numeric(c("0.001733102",  "0.681109185", 
                                                "0.002538071", "0.997461929", 
                                                "0.002958580", "0.623809524", 
                                                "2.170329670")),
        P.Stromal.Immune=as.numeric(c("0.09658928", "0.04585841", "0.67806841", 
                                      "0.32193159", "0.99704142", "0.25396825", 
                                      "0.23943162")))
    out <- calculate_proportions_of_cells_in_structure(
        spe_structure, cell_types_of_interest = c("Immune1","Immune3"),
        feature_colname = "Cell.Type")
    expect_equal(res, out, tolerance = 1e-4) 
    
    ## summarise the distances from cells to the tumour border
    res <- data.frame(
        Cell.Type = c("All_cell_types_of_interest", 
                      "All_cell_types_of_interest", 
                      "Immune1", "Immune1", "Immune3", "Immune3"),
        Area = c("Tumor_area", "Stroma", "Tumor_area", "Stroma", "Tumor_area", 
                 "Stroma"),
        Min_d=as.numeric(c("10.9322494547641", "10.0238703579346", NA, 
                           "84.2001833941197", "10.9322494547641", 
                           "10.0238703579346")),
        Max_d = as.numeric(c("192.409359800297", "971.56383420638", NA, 
                             "970.774932660564", "192.409359800297", 
                             "971.56383420638")),
        Mean_d = as.numeric(c("86.200421492396", "195.106365999404", NA,
                              "346.140958983386", "86.200421492396",
                              "102.79227480847")),
        Median_d=as.numeric(c("88.2329866304885", "101.951127102521", NA,
                              "301.015350157948", "88.2329866304885", 
                              "68.192180252124")),
        St.dev_d=as.numeric(c("45.27413945602", "194.685066632607",NA,
                              "187.042468626624", "45.27413945602", 
                              "131.327138494673")))
    out <- calculate_summary_distances_of_cells_to_borders(
        spe_structure,  cell_types_of_interest = c("Immune1","Immune3"),
        feature_colname = "Cell.Type")
    expect_equal(res, out) 
})

test_that("functions in spatial heterogeneity family work", {
    ## grid_metrics()
    res <- c(0.7793498, 0.9995910, 0.9612366, 0.9456603, 0.8812909,
             0.9915529, 0.7617327, 0.4586858, 0.9905577, 0.9798688)
    grid <- grid_metrics(defined_image, FUN = calculate_entropy, n_split = 5,
                         cell_types_of_interest=c("Tumour","Immune3"), 
                         feature_colname = "Cell.Type")
    # test if the data strcure is "RasterLayer"
    expect_s4_class(grid, "RasterLayer")   
    
    out <- grid@data@values[1:10]
    expect_equal(res, out, tolerance = 1e-4) # test if the results are the same
    
    ## calculate_percentage_of_grids()
    res <- 96
    out <- calculate_percentage_of_grids(grid, threshold = 0.75, above = TRUE)
    expect_equal(res, out, tolerance = 1e-4)
    
    ## calculate_spatial_autocorrelation
    res <- -0.1152657
    out <- calculate_spatial_autocorrelation(grid, metric = "globalmoran")
    expect_equal(res, out, tolerance = 1e-4)
})

test_that("functions in identify neighborhood family work", {
    
    ## identify_neighborhoods()
    neighborhoods <- identify_neighborhoods(
        image_no_markers, method = "hierarchical", min_neighborhood_size = 100,
        cell_types_of_interest = c("Immune", "Immune1", "Immune2"),
        radius = 50, feature_colname = "Cell.Type")
    
    # test if the data strcure is "SpatialExperiment"
    expect_s4_class(neighborhoods, "SpatialExperiment") 
    
    ## composition_of_neighborhoods()
    res <- data.frame(row.names = c(1L, 2L, 3L),
                      Cell.Type = c("Immune", "Immune1", "Immune2"), 
                      Neighborhood = rep("Cluster_1",3),
                      Number_of_cells = c(98, 143, 54),
                      Total_number_of_cells = rep(295,3),
                      Percentage = c(33.220339, 48.474576, 18.305085))
    
    neighborhoods_vis <- composition_of_neighborhoods(
        neighborhoods, feature_colname="Cell.Type")
    
    expect_equal(neighborhoods_vis[1:3, ], res, tolerance = 1e-4)
    
    ## plot_composition_heatmap()
    p <- plot_composition_heatmap(neighborhoods_vis, 
                                  feature_colname="Cell.Type")
    
    expect_is(p, "Heatmap")
})
