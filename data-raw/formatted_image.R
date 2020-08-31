## code to prepare `formatted_image` dataset

library(SPIAT)

raw_inform_data <- "S6_[49209,17530]_cell_seg_data.txt"
markers <- c("DAPI","CD3","PDL-1","CD4","CD8","AMACR")
locations <- c("Nucleus", "Cytoplasm", "Membrane", "Cytoplasm", "Cytoplasm", "Cytoplasm")

formatted_image <- format_image_to_sce(
                          format="INFORM",
                          image=raw_inform_data,
                          markers=markers,
                          locations= locations)

usethis::use_data(formatted_image, overwrite = TRUE)
