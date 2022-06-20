## code to prepare `simulated_images` dataset goes here
#####
library(spaSim)
library(SPIAT)

# simulate an image with three tumor clusters with immune rings.
set.seed(610)
bg <- simulate_mixing(idents = c("Tumour", "Immune1", "Immune2", 
                                          "Immune", "Others"),
                      props =c(0.03, 0.1, 0.05, 0.02, 0.8),
                      plot_image = FALSE)
# Tumour: Tumour_marker
# Immune1: Immune_marker1,Immune_marker2
# Immune2: Immune_marker1,Immune_marker3
# Immune: Immune_marker1,Immune_marker2,Immune_marker4
# Others: all other cells not included in above

set.seed(610)
a <- multiple_images_with_immune_rings(
    bg_sample = bg, cluster_size = 150, ring_width = 80,
    prop_ring_infiltration = 0.3, prop_infiltration = 0.1, cluster_loc_x = 50, 
    cluster_loc_y = -100, plot_image = FALSE
)
a <- a[[1]]
set.seed(610)
a <- multiple_images_with_immune_rings(
    bg_sample = a, ring_shape = 2, cluster_size = 120, ring_width = 80,
    prop_ring_infiltration = 0.3, prop_infiltration = 0.1, cluster_loc_x = -500, 
    cluster_loc_y = 500, plot_image = FALSE
)
a <- a[[1]]

plot_cell_categories(a, c("Tumour","Immune1", "Immune2", "Immune"), 
                     c("red","darkblue", "lightblue", "darkgreen"), "Cell.Type")

# simulate markers 
table(a$Cell.Type)
# Immune Immune1 Immune2  Others  Tumour 
# 630     338     178    2986     819 


n_cells <- sum(table(a$Cell.Type))
set.seed(610)
Tumour_marker <- rbeta(n_cells, 0.25, 0.78)
plot(density(Tumour_marker))
Tumour_marker_p <- Tumour_marker > 0.6
table(Tumour_marker_p)

set.seed(610)
Immune_marker1 <- rbeta(n_cells, 0.07, 0.5)
plot(density(Immune_marker1))
Immune_marker1_p <- Immune_marker1 > 0.2
table(Immune_marker1_p)

set.seed(610)
Immune_marker2 <- rbeta(n_cells, 0.08, 0.6)
plot(density(Immune_marker2))
Immune_marker2_p <- Immune_marker2 > 0.2
table(Immune_marker2_p)

set.seed(610)
Immune_marker3 <- rbeta(n_cells, 0.08, 3)
plot(density(Immune_marker3))
Immune_marker3_p <- Immune_marker3 > 0.2
table(Immune_marker3_p)

Immune_marker4 <- rbeta(n_cells, 1, 200)
plot(density(Immune_marker4))
Immune_marker4_p <- Immune_marker4 > 0.01
table(Immune_marker4_p)


# assign the marker intensities to cells
data <- get_colData(a)
data$sample_id <- NULL


# Tumour_marker
set.seed(610)
Tumour_marker_intensities <- sample(Tumour_marker[Tumour_marker_p], 819, 
                                    replace = TRUE)
data[data$Cell.Type == "Tumour", "Tumour_marker"] <- Tumour_marker_intensities
set.seed(610)
NoTumour_marker_intensities <- sample(Tumour_marker[!Tumour_marker_p], 
                                      n_cells-819, replace = TRUE)
data[data$Cell.Type != "Tumour", "Tumour_marker"] <- NoTumour_marker_intensities

#Immune_marker1
set.seed(610)
Immune_marker1_intensities <- sample(Immune_marker1[Immune_marker1_p], 1146, 
                                     replace = TRUE)
data[data$Cell.Type %in% c("Immune","Immune1","Immune2"), 
     "Immune_marker1"] <- Immune_marker1_intensities
set.seed(610)
NoImmune_marker1_intensities <- sample(Immune_marker1[!Immune_marker1_p], 
                                       n_cells-1146, replace = TRUE)
data[!data$Cell.Type %in% c("Immune","Immune1","Immune2"), 
     "Immune_marker1"] <- NoImmune_marker1_intensities

#Immune_marker2
set.seed(610)
Immune_marker2_intensities <- sample(Immune_marker2[Immune_marker2_p], 968, 
                                     replace = TRUE)
data[data$Cell.Type %in% c("Immune","Immune1"), 
     "Immune_marker2"] <- Immune_marker2_intensities
set.seed(610)
NoImmune_marker2_intensities <- sample(Immune_marker2[!Immune_marker2_p], 
                                       n_cells-968, replace = TRUE)
data[!data$Cell.Type %in% c("Immune","Immune1"), 
     "Immune_marker2"] <- NoImmune_marker2_intensities

#Immune_marker3
set.seed(610)
Immune_marker3_intensities <- sample(Immune_marker3[Immune_marker3_p], 178, 
                                     replace = TRUE)
data[data$Cell.Type == "Immune2", "Immune_marker3"] <- Immune_marker3_intensities
set.seed(610)
NoImmune_marker3_intensities <- sample(Immune_marker3[!Immune_marker3_p], 
                                       n_cells-178, replace = TRUE)
data[data$Cell.Type != "Immune2", "Immune_marker3"] <- NoImmune_marker3_intensities

#Immune_marker4
set.seed(610)
Immune_marker4_intensities <- sample(Immune_marker4[Immune_marker4_p], 630, 
                                     replace = TRUE)
data[data$Cell.Type == "Immune", "Immune_marker4"] <- Immune_marker4_intensities
set.seed(610)
NoImmune_marker4_intensities <- sample(Immune_marker4[!Immune_marker4_p], 
                                       n_cells-630, replace = TRUE)
data[data$Cell.Type != "Immune", "Immune_marker4"] <- NoImmune_marker4_intensities

# define the original phenotypes
data[data$Cell.Type == "Tumour", "Phenotype_ori"] <- "Tumour_marker"
data[data$Cell.Type == "Immune1", "Phenotype_ori"] <- "Immune_marker1,Immune_marker2"
data[data$Cell.Type == "Immune2", "Phenotype_ori"] <- "Immune_marker1,Immune_marker3"
data[data$Cell.Type == "Immune", "Phenotype_ori"] <- "Immune_marker1,Immune_marker2,Immune_marker4"
data[data$Cell.Type == "Others", "Phenotype_ori"] <- "OTHER"

# format to sp
intensity_matrix <- t(data[,c("Tumour_marker","Immune_marker1","Immune_marker2",
                              "Immune_marker3","Immune_marker4")])
colnames(intensity_matrix) <- data$Cell.ID
coord_x <- data$Cell.X.Position
coord_y <- data$Cell.Y.Position
phenotypes <- data$Phenotype_ori

simulated_image <- format_image_to_spe(format = "general",
                                       intensity_matrix = intensity_matrix,
                                       coord_x = coord_x,
                                       coord_y = coord_y,
                                       phenotypes = phenotypes)

#####
usethis::use_data(simulated_image, overwrite = TRUE)

