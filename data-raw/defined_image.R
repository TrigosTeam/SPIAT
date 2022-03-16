## code to prepare `defined_image` dataset goes here
#####
categories <- c("Tumour_marker", "Immune_marker1,Immune_marker2", 
               "Immune_marker1,Immune_marker3", 
               "Immune_marker1,Immune_marker2,Immune_marker4", "OTHER")
names <- c("Tumour", "Immune1", "Immune2", "Immune3", "Others")

defined_image <- define_celltypes(SPIAT::simulated_image, categories,
                                  category_colname = "Phenotype",
                                  names, new_colname = "Cell.Type")
#####
usethis::use_data(defined_image, overwrite = TRUE)
