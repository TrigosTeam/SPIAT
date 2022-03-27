#' Dimensionality reduction plot
#'
#' @description Generates the dimensionality reduction plots (UMAP or tSNE)
#'   based on marker intensities. Cells are grouped by the categories under the
#'   selected column. -- have you tried doing PCA on the matrix and then doing the UMAP/tSNE? Does it help? scRNAseq workflows do this.
#' @param sce_object SingleCellExperiment object in the form of the output of
#'   \code{\link{format_image_to_sce}}.
#' @param plot_type String. Choose from "UMAP" and "TSNE".
#' @param scale Boolean. Whether scale the marker intensities.
#' @param feature_colname String. Specify the column name to group the cells.
#'
#' @return A plot
#' @export
#' @examples 
#' dimensionality_reduction_plot(SPIAT::simulated_image, plot_type = "TSNE", 
#' feature_colname = "Phenotype")
dimensionality_reduction_plot <- function(sce_object, plot_type = "UMAP", 
                                          scale=TRUE, feature_colname){
    
    Cell_ID <- X_coord <- Y_coord <- Label <- NULL
    formatted_data <- get_colData(sce_object)
    
    intensity_matrix <- SummarizedExperiment::assay(sce_object)
    intensity_matrix_no_DAPI <- intensity_matrix[rownames(intensity_matrix) != "DAPI",]
    
    if(scale){
        intensity_matrix_no_DAPI_scaled <- scale(t(intensity_matrix_no_DAPI))
    }else{
        intensity_matrix_no_DAPI_scaled <- t(intensity_matrix_no_DAPI)
    }
    
    if(plot_type == "UMAP"){
        intensity_DR <- umap::umap(intensity_matrix_no_DAPI_scaled)
        intensity_DR_layout <- as.data.frame(intensity_DR$data)
        colnames(intensity_DR_layout) <- c("X_coord", "Y_coord")
        intensity_DR_layout$Label <- formatted_data[[feature_colname]][match(rownames(intensity_DR_layout), formatted_data$Cell.ID)]
        
    }else if (plot_type == "TSNE"){
        intensity_DR <- Rtsne::Rtsne(intensity_matrix_no_DAPI_scaled)
        intensity_DR_layout <- as.data.frame(intensity_DR$Y)
        colnames(intensity_DR_layout) <- c("X_coord", "Y_coord")
        rownames(intensity_DR_layout) <- rownames(intensity_matrix_no_DAPI_scaled)
        intensity_DR_layout$Label <- formatted_data[[feature_colname]][match(rownames(intensity_DR_layout), formatted_data$Cell.ID)]
        
    }else{
        print("Print select UMAP or TSNE as plot type")
    }
    
    intensity_DR_layout$Cell_ID <- rownames(intensity_DR_layout)
    
    g <- ggplot(intensity_DR_layout, aes(x=X_coord, y=Y_coord, label=Cell_ID))+
        geom_point(aes(colour=Label))+
        ggtitle(plot_type)+
        theme_bw()+
        theme(panel.grid.major = element_blank())
    plotly::ggplotly(g)
    
    return(g)
}
