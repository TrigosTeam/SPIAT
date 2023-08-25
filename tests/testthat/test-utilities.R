context("utilities")


test_that("format_image_to_spe() works as general format", {
    
    #Construct a marker intensity matrix (rows are markers, columns are cells)
    intensity_matrix <- matrix(c(14.557, 0.169, 1.655, 0.054, 17.588, 0.229,
                                 1.188, 2.074, 21.262, 4.206,  5.924, 0.021), 
                               nrow = 4, ncol = 3)
    # define marker names as rownames
    rownames(intensity_matrix) <- c("DAPI", "CD3", "CD4", "AMACR")
    # define cell IDs as colnames
    colnames(intensity_matrix) <- c("Cell_1", "Cell_2", "Cell_3")
    # Construct a dummy metadata (phenotypes, x/y coordinates)
    # the order of the elements in these vectors correspond to the cell order
    # in `intensity matrix`
    phenotypes <- c("OTHER",  "AMACR", "CD3,CD4")
    coord_x <- c(82, 171, 184)
    coord_y <- c(30, 22, 38)
    
    formatted_image <- format_image_to_spe(format = "general",
                                           intensity_matrix = intensity_matrix, 
                                           phenotypes = phenotypes, 
                                           coord_x = coord_x, coord_y = coord_y)
    
    expect_is(formatted_image, "SpatialExperiment")
    expect_equal(dim(formatted_image), c(4L, 3L))
})

test_that("format_image_to_spe() works for inForm format", {
    
    raw_inform_data <- system.file("extdata", "tiny_inform.txt.gz", 
                                   package = "SPIAT")
    markers <- c("DAPI","CD3","PD-L1","CD4","CD8","AMACR")
    intensity_columns_interest <- c(
        "Nucleus DAPI (DAPI) Mean (Normalized Counts, Total Weighting)",
        "Cytoplasm CD3 (Opal 520) Mean (Normalized Counts, Total Weighting)",
        "Membrane PD-L1 (Opal 540) Mean (Normalized Counts, Total Weighting)",
        "Cytoplasm CD4 (Opal 620) Mean (Normalized Counts, Total Weighting)",
        "Cytoplasm CD8 (Opal 650) Mean (Normalized Counts, Total Weighting)",
        "Cytoplasm AMACR (Opal 690) Mean (Normalized Counts, Total Weighting)"
    )
    spe <- format_image_to_spe(format = "inForm",
                               path=raw_inform_data,
                               markers=markers,
                               intensity_columns_interest=
                                   intensity_columns_interest)
    
    expect_is(spe, "SpatialExperiment")
    expect_equal(dim(spe), c(6L, 9L))
})

test_that("format_inform_to_spe() works", {
    
    raw_inform_data <- system.file("extdata", "tiny_inform.txt.gz", 
                                   package = "SPIAT")
    markers <- c("DAPI","CD3","PD-L1","CD4","CD8","AMACR")
    intensity_columns_interest <- c(
        "Nucleus DAPI (DAPI) Mean (Normalized Counts, Total Weighting)",
        "Cytoplasm CD3 (Opal 520) Mean (Normalized Counts, Total Weighting)",
        "Membrane PD-L1 (Opal 540) Mean (Normalized Counts, Total Weighting)",
        "Cytoplasm CD4 (Opal 620) Mean (Normalized Counts, Total Weighting)",
        "Cytoplasm CD8 (Opal 650) Mean (Normalized Counts, Total Weighting)",
        "Cytoplasm AMACR (Opal 690) Mean (Normalized Counts, Total Weighting)"
    )
    spe <- format_inform_to_spe(path=raw_inform_data, markers=markers,
                                intensity_columns_interest=
                                    intensity_columns_interest)
    
    expect_is(spe, "SpatialExperiment")
    expect_equal(dim(spe), c(6L, 9L))
})

test_that("format_halo_to_spe() works", {
    
    raw_halo_data <- system.file("extdata", "tiny_halo.csv.gz", package="SPIAT")
    markers <- c("DAPI", "CD3", "PDL-1", "CD4", "CD8", "AMACR")
    intensity_columns_interest <- 
        c("Dye 1 Nucleus Intensity", "Dye 2 Cytoplasm Intensity",
          "Dye 3 Membrane Intensity", "Dye 4 Cytoplasm Intensity", 
          "Dye 5 Cytoplasm Intensity", "Dye 6 Cytoplasm Intensity")
    dye_columns_interest <-
        c("Dye 1 Positive Nucleus","Dye 2 Positive Cytoplasm",
          "Dye 3 Positive Membrane", "Dye 4 Positive Cytoplasm",
          "Dye 5 Positive Cytoplasm", "Dye 6 Positive Cytoplasm")
    formatted_HALO <- format_halo_to_spe(
        path=raw_halo_data, markers=markers,
        intensity_columns_interest=intensity_columns_interest, 
        dye_columns_interest=dye_columns_interest)
    
    expect_is(formatted_HALO, "SpatialExperiment")
    expect_equal(dim(formatted_HALO), c(6L, 10L))
})

test_that("format_codex_to_spe() works", {
    
    path <- system.file("extdata", "tiny_codex.csv.gz", package = "SPIAT")
    path_to_codex_cell_phenotypes <- 
        system.file("extdata", "tiny_codex_phenotypes.txt.gz", 
                    package = "SPIAT")
    markers <- c("CD45", "Ly6C", "CD27", "CD5", "CD79b")
    formatted_codex <- format_codex_to_spe(path = path, markers = markers, 
                                           path_to_codex_cell_phenotypes = 
                                               path_to_codex_cell_phenotypes)
    
    expect_is(formatted_codex, "SpatialExperiment")
    expect_equal(dim(formatted_codex), c(8L, 9L))
})


test_that("format_cellprofiler_to_spe() works", {
    
    path <- system.file("extdata", "tiny_cellprofiler.txt.gz", 
                        package = "SPIAT")
    markers <- c("Marker1", "Marker2", "Marker3", "Marker4", "Marker5", 
                 "DAPI", "Marker6")
    intensity_columns_interest <- c("Intensity_MeanIntensity_Marker1_rs",
                                    "Intensity_MeanIntensity_Marker2_rs",
                                    "Intensity_MeanIntensity_Marker3_rs",
                                    "Intensity_MeanIntensity_Marker4_rs",
                                    "Intensity_MeanIntensity_Marker5_rs",
                                    "Intensity_MeanIntensity_DAPI_rs", 
                                    "Intensity_MeanIntensity_Marker6_rs")
    formatted_cellprofiler <- 
        format_cellprofiler_to_spe(path = path, markers = markers, 
                                   intensity_columns_interest = 
                                       intensity_columns_interest)
    
    expect_is(formatted_cellprofiler, "SpatialExperiment")
    expect_equal(dim(formatted_cellprofiler), c(7L, 9L))
})

test_that("format_spe_to_ppp() works", {
    
    ppp_object <- format_spe_to_ppp(defined_image, 
                                    feature_colname = "Cell.Type")
    
    expect_is(ppp_object, "ppp")
    expect_equal(ppp_object$x[1], 139.8, tolerance = 0.1)
    expect_equal(ppp_object$y[1], 86.7, tolerance = 0.1)
    expect_equal(ppp_object$marks[1], "Others")
    
})

test_that("select_phenotypes works", {
    
    phenotypes_keep <- c("Immune_marker1,Immune_marker2", "Tumour_marker")
    
    data_subset <- 
        select_celltypes(simulated_image, 
                         celltypes = c("Tumour_marker",
                                       "Immune_marker1,Immune_marker2"),
                         feature_colname = "Phenotype", keep=TRUE)
    phenotypes <- unique(colData(data_subset)$Phenotype)
    
    expect_equal(phenotypes, phenotypes_keep)
})

test_that("define_celltypes works", {
    
    celltypes_defined <- c("Others", "Immune1", "Tumour", "Immune3", "Immune2")
    
    defined_spe <- define_celltypes(
        simulated_image, 
        categories = c("Tumour_marker", "Immune_marker1,Immune_marker2", 
                       "Immune_marker1,Immune_marker3", 
                       "Immune_marker1,Immune_marker2,Immune_marker4", "OTHER"), 
        category_colname = "Phenotype",
        names = c("Tumour", "Immune1", "Immune2","Immune3", "Others"), 
        new_colname = "Cell.Type")
    celltypes <- unique(colData(defined_spe)$Cell.Type)
    
    expect_equal(celltypes, celltypes_defined)
})

test_that("image_splitter works", {
    
    res <- data.frame(
        row.names = c(1L, 2L, 3L, 4L),
        Cell.ID = c("Cell_1", "Cell_2", "Cell_3", "Cell_4"),
        Phenotype = c("OTHER", "OTHER", "OTHER", "OTHER"),
        Cell.X.Position = c(139.77484, 77.86721, 84.44626, 110.19857),
        Cell.Y.Position = c(86.704079, 80.096527, 19.238638, 5.656004))
    
    split_image <- image_splitter(simulated_image, number_of_splits=3, 
                                  plot = FALSE)
    out1 <- data.frame(colData(split_image[[1]])[1:4,])
    out2 <- data.frame(spatialCoords(split_image[[1]])[1:4,])
    out <- cbind(out1, out2)
    out$sample_id <- NULL
    rownames(out) <- c(1L, 2L, 3L, 4L)
    
    expect_is(split_image, "list")
    expect_is(split_image[[1]], "SpatialExperiment")
    expect_equal(length(split_image), 9)
    expect_equal(out, res)
})

test_that("predict_phenotypes works", {
    
    # test it works when reference_phenotypes = TRUE
    res <-    data.frame(
        row.names = c(1L, 2L),
        Cell.ID = c("Cell_1","Cell_2"),
        Tumour_marker_actual_phenotype = c(0,0),
        Immune_marker1_actual_phenotype = c(0,0),
        Immune_marker2_actual_phenotype = c(0,0),
        Immune_marker3_actual_phenotype = c(0,0),
        Immune_marker4_actual_phenotype = c(0,0),
        Immune_marker1_predicted_phenotype = c(0,0),
        Immune_marker2_predicted_phenotype = c(0,0),
        Immune_marker3_predicted_phenotype = c(0,0),
        Immune_marker4_predicted_phenotype = c(0,0),
        Tumour_marker_predicted_phenotype = c(1,0))
    
    predicted_result <- predict_phenotypes(
        spe_object = simulated_image, thresholds = NULL, 
        tumour_marker = "Tumour_marker",
        baseline_markers = c("Immune_marker1", "Immune_marker2", 
                             "Immune_marker3", "Immune_marker4"), 
        reference_phenotypes = TRUE)
    
    out <- predicted_result[1:2,c(1,11:20)]
    
    expect_is(predicted_result, "data.frame")
    expect_equal(dim(predicted_result), c(4951, 20))
    expect_equal(out, res)
    
    # test it works when reference_phenotypes = FALSE
    predicted_spe_image <- predict_phenotypes(
        spe_object = simulated_image, thresholds = NULL, 
        tumour_marker = "Tumour_marker",
        baseline_markers = c("Immune_marker1", "Immune_marker2", 
                             "Immune_marker3", "Immune_marker4"), 
        reference_phenotypes = FALSE)
    
    predicted_phenotypes <- unique(predicted_spe_image$Phenotype)
    
    expect_length(predicted_phenotypes, 27)
})
