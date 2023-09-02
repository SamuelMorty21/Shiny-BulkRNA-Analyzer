create_custom_heatmap <- function(data, col_data, genes_info, show_column_names, rect_gp = gpar(col = NA)) {
  heatmap_args <- list(
    data = data,
    name = "log2(TPM+1)",
    show_column_names = show_column_names,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    column_names_gp = gpar(fontface = "bold"),
    row_gap = unit(3, "mm"),
    column_names_rot = -60,
    col = col_data,
    row_split = factor(genes_info$Category, level = unique(genes_info$Category)),
    rect_gp = rect_gp
  )
  
  heatmap_object <- Heatmap(
    heatmap_args$data,
    name = heatmap_args$name,
    show_column_names = heatmap_args$show_column_names,
    cluster_rows = heatmap_args$cluster_rows,
    cluster_columns = heatmap_args$cluster_columns,
    show_row_names = heatmap_args$show_row_names,
    column_names_gp = heatmap_args$column_names_gp,
    row_gap = heatmap_args$row_gap,
    column_names_rot = heatmap_args$column_names_rot,
    col = heatmap_args$col,
    row_split = heatmap_args$row_split,
    rect_gp = heatmap_args$rect_gp
  )
  
  return(heatmap_object)
}