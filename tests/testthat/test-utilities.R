context("utilities")

test_that("format_image_to_sce() works", {
  
  load("testdata/test_sce.rda")
  
  
  raw_inform_data <- "testdata/test_spiat.txt"
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
    
    res <- c("AMACR", "CD3,CD4", "CD3,CD8", "PDL-1")
    
    phenotypes <- print_phenotypes(formatted_image)
  
    expect_equal(phenotypes, res)
})


test_that("select_phenotypes works", {
    
    phenotypes <- select_phenotypes(formatted_image, keep=TRUE,
                                         phenotypes = c("AMACR",
                                                        "CD3,CD4",
                                                        "CD3,CD8",
                                                        "CD3,CD8",
                                                        "PDL-1",
                                                        "AMACR,PDL-1"))
    
    expect_equal(phenotypes, formatted_image)
})


test_that("image_splitter works", {
    
    res <-    data.frame(row.names = c("Cell_2", "Cell_3", "Cell_4", "Cell_5"),
                         Phenotype = c("AMACR", "AMACR", "AMACR", "AMACR"
                         ), Cell.X.Position = c(171L, 184L, 201L, 219L), 
                         Cell.Y.Position = c(22L, 38L, 52L, 63L), 
                         Cell.Area = c(464L, 553L, 462L, 876L), 
                         Nucleus.Area = c(177L, 212L, 239L, 451L), 
                         Nucleus.Compactness = c(0.54, 0.51, 0.53, 0.53), 
                         Nucleus.Axis.Ratio = c(1.84, 1.92, 1.47, 1.19), 
                         Cell.Axis.Ratio = c(1, 1.21, 1.09, 1.34), 
                         Cell_type = c("AMACR", "AMACR", "AMACR", "AMACR"))
    
    split_image <- image_splitter(formatted_image, number_of_splits=3)
    
    expect_is(split_image, "list")
    expect_equal(length(split_image), 9)
    expect_equal(split_image[[1]][1:4, ], res)
})

