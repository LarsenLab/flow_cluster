


plot_clustering_heatmap_wrapper <- function(expr, expr01, 
                                            cell_clustering, color_clusters, cluster_merging = NULL) {
    
    # calculate median expression
    expr_median <- data.frame(expr, cell_clustering = cell_clustering) %>%
        group_by(cell_clustering) %>% 
        summarize_all(funs(median))
    expr01_median <- data.frame(expr01, cell_clustering = cell_clustering) %>%
        group_by(cell_clustering) %>% 
        summarize_all(funs(median))
    
    # calculate cluster frequency
    clustering_table <- as.numeric(table(cell_clustering))
    
    # this clustering is based on the markers that were used for the main clustering
    d <- dist(expr_median[, colnames(expr)], method = "euclidean")
    cluster_rows <- hclust(d, method = "average")
    
    expr_heat <- as.matrix(expr01_median[, colnames(expr01)])
    rownames(expr_heat) <- expr01_median$cell_clustering
    
    
    labels_row <- paste0(rownames(expr_heat))
    labels_col <- colnames(expr_heat)
    
    # row annotation for the heatmap
    annotation_row <- data.frame(cluster = factor(expr01_median$cell_clustering))
    rownames(annotation_row) <- rownames(expr_heat)
    
    color_clusters <- color_clusters[1:nlevels(annotation_row$cluster)]
    names(color_clusters) <- levels(annotation_row$cluster)
    annotation_colors <- list(cluster = color_clusters)
    annotation_legend <- FALSE
    
    
    # color for the heatmap
    color <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100)
    
    pheatmap(expr_heat, color = color, 
             cluster_cols = T, cluster_rows = cluster_rows, 
             labels_col = labels_col, labels_row = labels_row, 
             display_numbers = FALSE, number_color = "black", 
             fontsize = 8, fontsize_number = 4,
             annotation_row = annotation_row, annotation_colors = annotation_colors, 
             annotation_legend = annotation_legend, border_color = NA)
    
}


save_pheatmap_pdf <- function(x, filename, width = 8, height = 8) {
    stopifnot(!missing(x))
    stopifnot(!missing(filename))
    pdf(filename, width = width, height = height)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
}

