context("utilities")

test_that("format_image_to_sce() works for inForm data", {
  
  raw_inform_data <- system.file("extdata", "tiny_inform.txt.gz", package = "SPIAT")
  markers <- c("DAPI","CD3","PDL-1","CD4","CD8","AMACR")
  intensity_columns_interest <- c(
    "Nucleus DAPI (DAPI) Mean (Normalized Counts, Total Weighting)",
    "Cytoplasm CD3 (Opal 520) Mean (Normalized Counts, Total Weighting)", 
    "Membrane PDL-1 (Opal 540) Mean (Normalized Counts, Total Weighting)",
    "Cytoplasm CD4 (Opal 620) Mean (Normalized Counts, Total Weighting)",
    "Cytoplasm CD8 (Opal 650) Mean (Normalized Counts, Total Weighting)", 
    "Cytoplasm AMACR (Opal 690) Mean (Normalized Counts, Total Weighting)"
  )
  sce <- format_image_to_sce(
    format="INFORM",
    path=raw_inform_data,
    markers=markers,
    intensity_columns_interest=intensity_columns_interest,
    dye_columns_interest=NULL
  )
  
  expect_is(sce, "SingleCellExperiment")
  expect_equal(dim(sce), c(6L, 9L))
})

test_that("format_sce_to_ppp() wroks", {
  
  ppp_object <- format_sce_to_ppp(defined_image, feature_colname = "Cell.Type")
  
  expect_is(ppp_object, "ppp")
  expect_equal(ppp_object$x[1], 139.8, tolerance = 0.1)
  expect_equal(ppp_object$y[1], 86.7, tolerance = 0.1)
  expect_equal(ppp_object$marks[1], "Others")
  
})

test_that("format_colData_to_sce() wroks", {
  
  df <- data.frame(Cell.ID = c("Cell_1", "Cell_2"), Cell.X.Positions = c(2,5), 
                   Cell.Y.Positions = c(3.3, 8), Phenotypes = c("CD3", "CD3,CD8"))
  sce <- format_colData_to_sce(df)
  
  out <- data.frame(SummarizedExperiment::colData(sce))
  
  expect_s4_class(sce, "SingleCellExperiment")
  expect_equal(dim(sce)[1], 1)
  expect_equal(dim(out)[2], 5)
})

test_that("print_feature works", {
    
    res <- c("OTHER","Immune_marker1,Immune_marker2", "Tumour_marker" , 
             "Immune_marker1,Immune_marker2,Immune_marker4", "Immune_marker1,Immune_marker3"   )
    
    phenotypes <- print_feature(simulated_image, "Phenotype")
  
    expect_equal(phenotypes, res)
})

test_that("select_phenotypes works", {
    
    phenotypes_keep <- c("Immune_marker1,Immune_marker2", "Tumour_marker")
    
    data_subset <- select_celltypes(simulated_image,
                                    celltypes = c("Tumour_marker","Immune_marker1,Immune_marker2"),
                                    feature_colname = "Phenotype", keep=TRUE)
    phenotypes <- unique(colData(data_subset)$Phenotype)
    
    expect_equal(phenotypes, phenotypes_keep)
})

test_that("define_celltypes works", {
  
  celltypes_defined <- c("Others", "Immune1", "Tumour", "Immune3", "Immune2")
  
  defined_sce <- define_celltypes(simulated_image, 
                                  categories = c("Tumour_marker", "Immune_marker1,Immune_marker2", "Immune_marker1,Immune_marker3", 
                                                 "Immune_marker1,Immune_marker2,Immune_marker4", "OTHER"), category_colname = "Phenotype",
                                  names = c("Tumour", "Immune1", "Immune2","Immune3", "Others"), new_colname = "Cell.Type")
  celltypes <- unique(colData(defined_sce)$Cell.Type)
  
  expect_equal(celltypes, celltypes_defined)
})

test_that("image_splitter works", {

   res <-    data.frame(
     row.names = c("Cell_1", "Cell_2", "Cell_3", "Cell_4"),
     Phenotype = c("OTHER", "OTHER", "OTHER", "OTHER"),
     Cell.X.Position = c(139.77484, 77.86721, 84.44626, 110.19857),
     Cell.Y.Position = c(86.704079, 80.096527, 19.238638, 5.656004))

   split_image <- image_splitter(simulated_image, number_of_splits=3, plot = FALSE)
   
   out <- split_image[[1]][1:4, ]

   expect_is(split_image, "list")
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
  
  predicted_result <- predict_phenotypes(sce_object = simulated_image, thresholds = NULL,
                                         tumour_marker = "Tumour_marker",
                                         baseline_markers = c("Immune_marker1", "Immune_marker2", 
                                                              "Immune_marker3", "Immune_marker4"), 
                                         reference_phenotypes = TRUE)
  
  out <- predicted_result[1:2,c(1,10:19)]
  
  expect_is(predicted_result, "data.frame")
  expect_equal(dim(predicted_result), c(4951, 19))
  expect_equal(out, res)
  
  # test it works when reference_phenotypes = FALSE
  predicted_sce_image <- predict_phenotypes(sce_object = simulated_image, thresholds = NULL,
                                            tumour_marker = "Tumour_marker",
                                            baseline_markers = c("Immune_marker1", "Immune_marker2", 
                                                                 "Immune_marker3", "Immune_marker4"), 
                                            reference_phenotypes = FALSE)
  
  predicted_phenotypes <- unique(predicted_sce_image$Phenotype)
  
  expect_length(predicted_phenotypes, 27)
})
