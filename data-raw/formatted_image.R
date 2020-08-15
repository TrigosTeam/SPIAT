## code to prepare `formatted_image` dataset

library(SPIAT)

raw_inform_data <- "S6_[49209,17530]_cell_seg_data.txt"

markers <- c("DAPI","CD3","PDL1","FOXP3","CD4","CD8","AMACR")

intensity_columns_interest <- c("Nucleus DAPI (DAPI) Mean (Normalized Counts, Total Weighting)","Cytoplasm CD3 (Opal 520) Mean (Normalized Counts, Total Weighting)", 
                                "Membrane PDL-1 (Opal 540) Mean (Normalized Counts, Total Weighting)", "Cytoplasm FOXP3 (Opal 570) Mean (Normalized Counts, Total Weighting)", 
                                "Cytoplasm CD4 (Opal 620) Mean (Normalized Counts, Total Weighting)", "Cytoplasm CD8 (Opal 650) Mean (Normalized Counts, Total Weighting)", 
                                "Cytoplasm AMACR (Opal 690) Mean (Normalized Counts, Total Weighting)")

#Formats an INFORM or HALO image into a singlecellexperiment class
#where the count assay stores the expression level of every marker (rows) for
#every cell (columns), and cell phenotype, x and y coordinates are stored
#under colData
formatted_image <- format_image_to_sce(format="INFORM",
                                       image=raw_inform_data,
                                       markers=markers,
                                       dye_columns_interest=NULL,
                                       intensity_columns_interest=intensity_columns_interest)


usethis::use_data(formatted_image, overwrite = TRUE)
