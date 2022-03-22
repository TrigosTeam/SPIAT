context("families")

test_that("functions in cross K family work", {
    ## calculate_cross_functions()
    res <- data.frame(row.names = c(10L, 11L, 12L),
                      r = c(9.09090909090909171653, 10.10101010101010210462, 11.11111111111111249272), 
                      theo = c(259.63575649502428177584, 320.53797098151147793033,387.85094488762888431665),
                      border = c(0.00000000000000000000,7.75578194807568443991,77.55781948075684795185),
                      stringsAsFactors = FALSE)
    df_cross <- calculate_cross_functions(defined_image, method = "Kcross",
                    cell_types_of_interest = c("Tumour","Immune3"),
                    feature_colname ="Cell.Type", dist = 100, plot_results = FALSE)
     # test if the data strcure is "fv"
    expect_s3_class(df_cross, "fv")   
    
    out <- data.frame(df_cross[10:12,])
    expect_equal(res, out) # test if the results are the same
    
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
    sce_border <- identify_bordering_cells(defined_image, reference_cell = "Tumour",
                                           feature_colname = "Cell.Type", n_to_exclude = 10)
    # test if the data strcure is "SingleCellExperiment"
    expect_s4_class(sce_border, "SingleCellExperiment") 
    
    out <- sce_border$Region[3776:3779]
    expect_equal(res, out) # test if the results are the same
    
    ## R_BT()
    # test if R_BT returns a number
    n <- R_BT(defined_image, cell_type_of_interest = "Tumour", "Cell.Type")
    
    expect_is(n, "numeric")
    
    ## calculate_distance_to_tumour_margin()
    res <- c(83.54843614991348488275, 63.04062766766428183018, 
             0.00000000000000000000, 29.59984593965495847101)
    res1 <- c("Non-border", "Non-border", "Border", "Non-border")
    sce_dist <- calculate_distance_to_tumour_margin(sce_border)
    # test if the data strcure is "SingleCellExperiment"
    expect_s4_class(sce_dist, "SingleCellExperiment") 
    
    out <- sce_dist$Distance.To.Border[3776:3779]
    out1 <- sce_dist$region2[3776:3779]
    # test if the results are the same
    expect_equal(res, out) 
    expect_equal(res1, out1)
    
    ## define_structure()
    res <- c("Internal.margin", "Internal.margin", "Border", "Internal.margin")
    sce_structure <- define_structure(sce_dist, names_of_immune_cells = c("Immune1","Immune2","Immune3"),
                                      feature_colname = "Cell.Type", n_margin_layers = 5)
    
    # test if the data strcure is "SingleCellExperiment"
    expect_s4_class(sce_structure, "SingleCellExperiment") 
    
    out <- sce_structure$Structure[3776:3779]
    expect_equal(res, out)  # test if the results are the same
    
    ## summarise the cell proportions in each tumour structure
    res <- data.frame(Cell.Type = c("Immune1", "Immune3", "Immune1", "Immune3", "Immune1", "Immune3", "All_cells_of_interest"),
                      Relative_to = c("All_cells_in_the_structure", "All_cells_in_the_structure", "All_cells_of_interest_in_the_structure", "All_cells_of_interest_in_the_structure", "The_same_cell_type_in_the_whole_image", "The_same_cell_type_in_the_whole_image", "All_cells_in_the_structure"),
                      P.Infiltrated.Immune=as.numeric(c("0", "0.143859649122807", "0", "1", "0", "0.0650793650793651", "0.143859649122807")),
                      P.Internal.Margin.Immune = as.numeric(c("0", "0.0878048780487805", "0", "1", "0", "0.0571428571428571", "0.0878048780487805")),
                      P.External.Margin.Immune = as.numeric(c("0.00549450549450549", "2.15934065934066", "0.00253807106598985", "0.99746192893401", "0.0029585798816568", "0.623809523809524", "2.17032967032967")),
                      P.Stromal.Immune=as.numeric(c("0.119715808170515", "0.0568383658969805", "0.678068410462777", "0.321931589537223", "0.997041420118343", "0.253968253968254", "0.23943161634103")))
    out <- calculate_proportions_of_cells_in_structure(sce_structure, 
                                                cell_types_of_interest = c("Immune1","Immune3"),
                                                feature_colname = "Cell.Type")
    expect_equal(res, out) 
    
    ## summarise the distances from cells to the tumour border
    res <- data.frame(Cell.Type = c("All_cell_types_of_interest", "All_cell_types_of_interest", "Immune1", "Immune1", "Immune3", "Immune3"),
                      Area = c("Tumor_area", "Stroma", "Tumor_area", "Stroma", "Tumor_area", "Stroma"),
                      Min_d=as.numeric(c("10.9322494547641", "10.0238703579346", Inf, "84.2001833941197", "10.9322494547641", "10.0238703579346")),
                      Max_d = as.numeric(c("192.409359800297", "971.56383420638", -Inf, "970.774932660564", "192.409359800297", "971.56383420638")),
                      Mean_d = as.numeric(c("86.200421492396", "195.106365999404", NaN, "346.140958983386", "86.200421492396", "102.79227480847")),
                      Median_d=as.numeric(c("88.2329866304885", "101.951127102521", NA, "301.015350157948", "88.2329866304885", "68.192180252124")),
                      St.dev_d=as.numeric(c("45.27413945602", "194.685066632607",NA, "187.042468626624", "45.27413945602", "131.327138494673")))
    out <- calculate_summary_distances_of_cells_to_borders(sce_structure, 
                                                           cell_types_of_interest = c("Immune1","Immune3"),
                                                           feature_colname = "Cell.Type")
    expect_equal(res, out) 
})

test_that("functions in spatial heterogeneity family work", {
    ## grid_metrics()
    res <- c(0.7793498, 0.9995910, 0.9612366, 0.9456603, 0.8812909,
             0.9915529, 0.7617327, 0.4586858, 0.9905577, 0.9798688)
    grid <- grid_metrics(defined_image, FUN = calculate_entropy, n_split = 5,
                         cell_types_of_interest=c("Tumour","Immune3"), feature_colname = "Cell.Type")
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
    res <- data.frame(row.names = c("Cell_15", "Cell_22", "Cell_23", "Cell_25"),
                      Cell.ID = c("Cell_15", "Cell_22", "Cell_23", "Cell_25"), 
                      Cell.X.Position = c(109.670273, 142.033630, 46.017590, 153.227951), 
                      Cell.Y.Position = c(17.129564, 32.068722, 4.760143, 128.299149),
                      Cell.Type = c("Immune1", "Immune", "Immune2", "Immune1"), 
                      Cell.Size = c(10.913455, 9.853069, 9.315296, 10.547904),
                      Cluster = rep("Free_cell",4))
    
    neighborhoods <- identify_neighborhoods(image_no_markers, method = "hierarchical",
                                            min_neighborhood_size = 100, 
                                            cell_types_of_interest = c("Immune", "Immune1", "Immune2"),
                                            radius = 50, feature_colname = "Cell.Type")
    
    expect_equal(neighborhoods[1:4, ], res, tolerance = 1e-1)
    
    ## composition_of_neighborhoods()
    res <- data.frame(row.names = c(1L, 2L, 3L),
                      Cell.Type = c("Immune", "Immune1", "Immune2"), 
                      Cluster = rep("Cluster_1",3),
                      Number_of_cells = c(98, 143, 54),
                      Total_number_of_cells = rep(295,3),
                      Percentage = c(33.220339, 48.474576, 18.305085))
    
    neighborhoods_vis <- composition_of_neighborhoods(neighborhoods, feature_colname="Cell.Type")
    
    expect_equal(neighborhoods_vis[1:3, ], res, tolerance = 1e-4)
    
    ## plot_composition_heatmap()
    p <- plot_composition_heatmap(neighborhoods_vis, feature_colname="Cell.Type")
    
    expect_is(p, "Heatmap")
})
