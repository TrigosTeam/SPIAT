context("utilities")

test_that("format_image_to_sce() works", {
  
  load("testdata/test_sce.rda")
  
  
  raw_inform_data <- raw_inform_data <- system.file("extdata", "tiny_inform.txt.gz", package = "SPIAT")
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
    image=raw_inform_data,
    markers=markers,
    intensity_columns_interest=intensity_columns_interest,
    dye_columns_interest=NULL
  )
  
  expect_equal(sce, test_sce)
})


test_that("print_phenotypes works", {
    
    res <- c("OTHER", "AMACR", "CD3,CD4", "CD3,CD8", "CD3", "PDL-1")
    
    phenotypes <- print_phenotypes(formatted_image)
  
    expect_equal(phenotypes, res)
})


test_that("select_phenotypes works", {
    
    phenotypes <- select_phenotypes(formatted_image, keep=TRUE,
                                         phenotypes = c("OTHER", "AMACR",
                                                        "CD3,CD4",
                                                        "CD3,CD8",
                                                        "CD3",
                                                        "PDL-1"))
    
    expect_equal(phenotypes, formatted_image)
})


test_that("image_splitter works", {
    
    res <-    data.frame(
      row.names = c("Cell_1", "Cell_2", "Cell_3", "Cell_4"),
      Phenotype = c("OTHER", "AMACR", "AMACR", "AMACR"), Cell.X.Position = c(82, 171, 184, 201), 
      Cell.Y.Position = c(30, 22, 38, 52), 
      Cell.Area = c(477, 464, 553, 462), 
      Nucleus.Area = c(160,177, 212, 239), 
      Nucleus.Compactness = c(0.52, 0.54, 0.51, 0.53), 
      Nucleus.Axis.Ratio = c(2.05, 1.84, 1.92, 1.47), 
      Cell.Axis.Ratio = c(1.48, 1, 1.21, 1.09), 
      Cell_type = c("OTHER", "AMACR", "AMACR", "AMACR"))
    
    split_image <- image_splitter(formatted_image, number_of_splits=3)
    
    expect_is(split_image, "list")
    expect_equal(length(split_image), 9)
    expect_equal(split_image[[1]][1:4, ], res)
})

