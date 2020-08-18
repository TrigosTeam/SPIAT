context("utilities")

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
