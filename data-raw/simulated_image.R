## code to prepare `DATASET` dataset goes here
#####
library(spaSim)
library(SPIAT)

# simulate an image with three tumor clusters with immune rings.
set.seed(610)
bg <- simulate_mixing(names_of_mixing = c("Tumour", "Immune1", "Immune2", 
                                          "Immune", "Others"),
                      mixing_degree =c(0.1, 0.1, 0.05, 0.02, 0.73),
                      plot.image = FALSE)
# Tumour: AMACR
# Immune1: CD3,CD4
# Immune2: CD3,CD8
# Immune: CD3,CD4,FOXP3
# Others: all other cells not included in above

set.seed(610)
a <- multiple_images_with_immune_rings(
    background_sample = bg, cluster_size = 150, ring_width = 80,
    ring_infiltration = 0.3, infiltration = 0.2, cluster_loc_x = 50, 
    cluster_loc_y = -100, plot.image = FALSE
)
a <- a[[1]]
set.seed(610)
a <- multiple_images_with_immune_rings(
    background_sample = a, ring_shape = 2, cluster_size = 120, ring_width = 80,
    ring_infiltration = 0.3, infiltration = 0.3, cluster_loc_x = -500, 
    cluster_loc_y = 500, plot.image = FALSE
)
a <- a[[1]]

a$Cell.Type <- a$Phenotype
a$Phenotype <- NULL
plot_cell_categories(a, c("Tumour","Immune1", "Immune2", "Immune"), 
                     c("red","darkblue", "lightblue", "darkgreen"), "Cell.Type")

# simulate markers 
table(a$Cell.Type)
#Immune Immune1 Immune2  Others  Tumour 
#757     359     159    2754     922 


n_cells <- sum(table(a$Cell.Type))
set.seed(610)
AMACR <- rbeta(n_cells, 0.25, 0.78)
plot(density(AMACR))
AMACR_p <- AMACR > 0.55
table(AMACR_p)

set.seed(610)
CD3 <- rbeta(n_cells, 0.7, 25)
plot(density(CD3))
CD3_p <- CD3 > 0.05
table(CD3_p)

set.seed(610)
CD4 <- rbeta(n_cells, 1, 20)
plot(density(CD4))
CD4_p <- CD4 > 0.07
table(CD4_p)

set.seed(610)
CD8 <- rbeta(n_cells, 1, 300)
plot(density(CD8))
CD8_p <- CD8 > 0.01
table(CD8_p)

FOXP3 <- rbeta(n_cells, 1, 200)
plot(density(FOXP3))
FOXP3_p <- FOXP3 > 0.01
table(FOXP3_p)


# assign the marker intensities to cells
data <- data.frame(colData(a))
data$lab <- NULL
data$pseudo <- NULL


# AMACR
set.seed(610)
AMACR_intensities <- sample(AMACR[AMACR_p], 922, replace = TRUE)
data[data$Cell.Type == "Tumour", "AMACR"] <- AMACR_intensities
set.seed(610)
NoAMACR_intensities <- sample(AMACR[!AMACR_p], n_cells-922, replace = TRUE)
data[data$Cell.Type != "Tumour", "AMACR"] <- NoAMACR_intensities

#CD3
set.seed(610)
CD3_intensities <- sample(CD3[CD3_p], 1275, replace = TRUE)
data[data$Cell.Type %in% c("Immune","Immune1","Immune2"), "CD3"] <- CD3_intensities
set.seed(610)
NoCD3_intensities <- sample(CD3[!CD3_p], n_cells-1275, replace = TRUE)
data[!data$Cell.Type %in% c("Immune","Immune1","Immune2"), "CD3"] <- NoCD3_intensities

#CD4
set.seed(610)
CD4_intensities <- sample(CD4[CD4_p], 1116, replace = TRUE)
data[data$Cell.Type %in% c("Immune","Immune1"), "CD4"] <- CD4_intensities
set.seed(610)
NoCD4_intensities <- sample(CD4[!CD4_p], n_cells-1116, replace = TRUE)
data[!data$Cell.Type %in% c("Immune","Immune1"), "CD4"] <- NoCD4_intensities

#CD8
set.seed(610)
CD8_intensities <- sample(CD8[CD8_p], 159, replace = TRUE)
data[data$Cell.Type == "Immune2", "CD8"] <- CD8_intensities
set.seed(610)
NoCD8_intensities <- sample(CD8[!CD8_p], n_cells-159, replace = TRUE)
data[data$Cell.Type != "Immune2", "CD8"] <- NoCD8_intensities

#FOXP3
set.seed(610)
FOXP3_intensities <- sample(FOXP3[FOXP3_p], 757, replace = TRUE)
data[data$Cell.Type == "Immune", "FOXP3"] <- FOXP3_intensities
set.seed(610)
NoFOXP3_intensities <- sample(FOXP3[!FOXP3_p], n_cells-757, replace = TRUE)
data[data$Cell.Type != "Immune", "FOXP3"] <- NoFOXP3_intensities

# define the orginial phenotypes
data[data$Cell.Type == "Tumour", "Phenotype_ori"] <- "AMACR"
data[data$Cell.Type == "Immune1", "Phenotype_ori"] <- "CD3,CD4"
data[data$Cell.Type == "Immune2", "Phenotype_ori"] <- "CD3,CD8"
data[data$Cell.Type == "Immune", "Phenotype_ori"] <- "CD3,CD4,FOXP3"
data[data$Cell.Type == "Others", "Phenotype_ori"] <- "OTHER"


# format to sce
intensity_matrix <- t(data[,c("AMACR","CD3","CD4","CD8","FOXP3")])
coord_x <- data$Cell.X.Position
coord_y <- data$Cell.Y.Position

simulated_image <- format_image_to_sce(format = "general",
                                       intensity_matrix = intensity_matrix,
                                       coord_x = coord_x,
                                       coord_y = coord_y,
                                       phenotypes = data$Phenotype_ori)

#####
usethis::use_data(simulated_image, overwrite = TRUE)
