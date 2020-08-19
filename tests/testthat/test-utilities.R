context("utilities")

test_that("format_image_to_sce() works", {
  
  load("testdata/test_sce.rda")
  
  
  raw_inform_data <- "testdata/test_spiat.txt"
  markers <- c("DAPI","CD3","PDL1","CD4","CD8","AMACR")
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
    
    res <- c("AMACR", "CD3,CD4", "CD3,CD8", "PDL1")
    
    phenotypes <- print_phenotypes(formatted_image)
  
    expect_equal(phenotypes, res)
})


test_that("select_phenotypes works", {
    
    phenotypes <- select_phenotypes(formatted_image, keep=TRUE,
                                         phenotypes = c("AMACR",
                                                        "CD3,CD4",
                                                        "CD3,CD8",
                                                        "CD3,CD8",
                                                        "PDL1",
                                                        "AMACR,PDL1"))
    
    expect_equal(phenotypes, formatted_image)
})


test_that("image_splitter works", {
    
    res <-    data.frame(row.names = c("Cell_2", "Cell_3", "Cell_4", "Cell_5"),
                         Phenotype = rep("AMACR", 4),
                         Cell.X.Position = c(171, 184, 201, 219),
                         Cell.Y.Position = c(22, 38, 52, 63),
                         Cell_type = rep("AMACR", 4))
    
    split_image <- image_splitter(formatted_image, number_of_splits=3)
    
    expect_is(split_image, "list")
    expect_equal(length(split_image), 9)
    expect_equal(split_image[[1]][1:4, ], res)
})


test_that("gg_color_hue() works", {
  
  res <- c("#F8766D", "#E18A00", "#BE9C00", "#8CAB00", "#24B700", 
           "#00BE70", "#00C1AB", "#00BBDA", "#00ACFC", "#8B93FF", 
           "#D575FE", "#F962DD", "#FF65AC")
  
  hex_colours <- gg_color_hue(13)
  
  expect_equal(hex_colours, res)
})


test_that("colhex2col() works", {
  
  res <- c("salmon", "orange3", "gold3", "yellow4", "green3", "springgreen3", 
           "turquoise3", "turquoise3", "deepskyblue2", "lightslateblue", 
           "mediumorchid1", "orchid2", "hotpink")
  
  hex_colours <- gg_color_hue(13)
  colour_vector <- sapply(hex_colours, colhex2col, USE.NAMES = FALSE)
  
  expect_equal(colour_vector, res)
})
