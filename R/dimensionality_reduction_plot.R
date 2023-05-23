#' Dimensionality reduction plot
#'
#' @description Generates the dimensionality reduction plots (UMAP or tSNE)
#'   based on marker intensities. Cells are grouped by the categories under the
#'   selected column.
#' @param spe_object SpatialExperiment object in the form of the output of
#'   \code{\link{format_image_to_spe}}.
#' @param plot_type String. Choose from "UMAP" and "TSNE".
#' @param scale Boolean. Whether scale the marker intensities.
#' @param perplexity Numeric. Perplexity parameter of the Rtsne function (should
#'   be positive and no bigger than 3 * perplexity < n - 1, where n is the
#'   number of cells).
#' @param feature_colname String. Specify the column name to group the cells.
#'
#' @return A plot
#' @export
#' @examples
#' dimensionality_reduction_plot(SPIAT::simulated_image, plot_type = "TSNE",
#' feature_colname = "Phenotype")
dimensionality_reduction_plot <- function(spe_object, plot_type = "UMAP", 
                                          scale=TRUE, 
                                          perplexity = 30, feature_colname){
    requireNamespace("plotly", quietly = TRUE)
    Cell_ID <- dim_X <- dim_Y <- Label <- NULL
    formatted_data <- get_colData(spe_object)
    
    intensity_matrix <- SummarizedExperiment::assay(spe_object)
    intensity_matrix_no_DAPI <- 
        intensity_matrix[rownames(intensity_matrix) != "DAPI",]
    
    if(scale){
        intensity_matrix_no_DAPI_scaled <- scale(t(intensity_matrix_no_DAPI))
    }else{
        intensity_matrix_no_DAPI_scaled <- t(intensity_matrix_no_DAPI)
    }
    
    if(plot_type == "UMAP"){
        if (dim(formatted_data)[1] <= 14){
            stop("The image should contain at least 15 cells to plot UMAP.")
        }
        requireNamespace("umap", quietly = TRUE)
        intensity_DR <- umap::umap(intensity_matrix_no_DAPI_scaled)
        intensity_DR_layout <- as.data.frame(intensity_DR$data)
        colnames(intensity_DR_layout) <- c("dim_X", "dim_Y")
        intensity_DR_layout$Label <- 
            formatted_data[[feature_colname]][
                match(rownames(intensity_DR_layout), formatted_data$Cell.ID)]
        
    }else if (plot_type == "TSNE"){
        if (dim(formatted_data)[1] <= 4){
            stop("The image should contain at least 5 cells to plot TSNE plot.")
        }
        if (3 * perplexity > nrow(formatted_data)- 1){
            perplexity <- floor((nrow(formatted_data)-1)/3)
            warning("The perplexity for tsne has been changed to ", perplexity,
                    " due to the small cell count in the image.")
        }
        requireNamespace("Rtsne", quietly = TRUE)
        intensity_DR <- Rtsne::Rtsne(intensity_matrix_no_DAPI_scaled,
                                     perplexity = perplexity)
        intensity_DR_layout <- as.data.frame(intensity_DR$Y)
        colnames(intensity_DR_layout) <- c("dim_X", "dim_Y")
        rownames(intensity_DR_layout) <- 
            rownames(intensity_matrix_no_DAPI_scaled)
        intensity_DR_layout$Label <- 
            formatted_data[[feature_colname]][
                match(rownames(intensity_DR_layout), formatted_data$Cell.ID)]
        
    }else{
        methods::show("Print select UMAP or TSNE as plot type")
    }
    
    intensity_DR_layout$Cell_ID <- rownames(intensity_DR_layout)
    
    g <- ggplot(intensity_DR_layout, aes(x=dim_X, y=dim_Y, label=Cell_ID))+
        geom_point(aes(colour=Label))+
        ggtitle(plot_type)+
        theme_bw()+
        theme(panel.grid.major = element_blank())
    plotly::ggplotly(g)
    
    return(g)
}
