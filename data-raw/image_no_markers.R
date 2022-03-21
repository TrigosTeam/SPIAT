## code to prepare `image_no_markers` dataset goes here
#####
library(spaSim)
library(SPIAT)

# simulate an image with three immune clusters
set.seed(610)
bg <- simulate_mixing(names_of_mixing = c("Tumour", "Immune1", "Immune2", 
                                          "Immune", "Others"),
                      mixing_degree =c(0.03, 0.1, 0.05, 0.02, 0.8),
                      plot.image = FALSE)

set.seed(610)
a <- simulate_clusters(background_sample = bg,
                       n_clusters = 4, properties_of_clusters = list(
                           C1 = list(name_of_cluster_cell = "Immune1", size = 900,
                                     shape = "Irregular", centre_loc = data.frame(x = 200, y = 200), 
                                     infiltration_types = c("Immune","Immune2", "Others"), 
                                     infiltration_proportions = c(0.3, 0.2, 0.1)),
                           C2 = list(name_of_cluster_cell = "Immune2", size = 1000,
                                     shape = "Irregular", centre_loc = data.frame(x = 1200, y = 400), 
                                     infiltration_types = c("Immune", "Others"), 
                                     infiltration_proportions = c(0.2, 0.05)),
                           C3 = list(name_of_cluster_cell = "Immune", size = 800, 
                                     shape = "Irregular", centre_loc = data.frame(x = 500, y = 1200), 
                                     infiltration_types = c("Immune1", "Immune2", "Others"), 
                                     infiltration_proportions = c(0.2, 0.3, 0.05)),
                           C4 = list(name_of_cluster_cell = "Immune", size = 280, 
                                     shape = "Oval", centre_loc = data.frame(x = 800, y = 1200), 
                                     infiltration_types = c("Immune1", "Immune2", "Others"), 
                                     infiltration_proportions = c(0.2, 0.3, 0.05))),
                       plot.image = TRUE)


a$Cell.Type <- a$Phenotype
a$Phenotype <- NULL
plot_cell_categories(a, c("Tumour","Immune1", "Immune2", "Immune"), 
                     c("red","darkblue", "darkgreen", "brown"), "Cell.Type")


table(a$Cell.Type)
# Immune Immune1 Immune2  Others  Tumour 
# 499     600     674    3059     119 
sum(table(a$Cell.Type))


# simulate cell sizes
set.seed(610)
Tumour_size <- rnorm(119, 15, 2)
plot(density(Tumour_size))
a[a$Cell.Type == "Tumour", "Cell.Size"] <- Tumour_size

set.seed(610)
Immune_size <- rnorm(1773, 9, 3)
plot(density(Immune_size))
a[a$Cell.Type %in% c("Immune","Immune1","Immune2"), "Cell.Size"] <- Immune_size

set.seed(610)
Other_size <- rnorm(3059, 12, 6)
plot(density(Other_size))
a[a$Cell.Type == "Others", "Cell.Size"] <- Other_size

# format to sce
intensity_matrix <- NULL
coord_x <- a$Cell.X.Position
coord_y <- a$Cell.Y.Position

image_no_markers <- format_colData_to_sce(a)
plot_cell_categories(image_no_markers, c("Tumour","Immune1", "Immune2", "Immune"), 
                     c("red","darkblue", "darkgreen", "brown"), "Cell.Type")
#####
usethis::use_data(image_no_markers, overwrite = TRUE)
