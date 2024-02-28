###################################
### PAGA connectivity heat map ####
###################################

paga_conn_hm <- function(mat, column_title="PAGA connectivity", name="Conect.", diff=FALSE, color_min=rev(brewer.pal(11,"RdBu"))[1], color_max=rev(brewer.pal(11,"RdBu"))[11]) {
    
    # Set color gradient 
    if(diff) {
        
        col=colorRamp2(c(-1, 0, 1), c(color_min, "#ffffff", color_max))
        
    } else {
        
        col=colorRamp2(c(0, 1), c("#ffffff", color_max))
    }
    
    hm <- ComplexHeatmap::Heatmap(
    
        matrix=mat,
        column_title=column_title,
        name=name, 
        column_title_gp=gpar(fontsize=10, fontface="bold"), 
        column_names_gp =grid::gpar(fontsize=6), 
        row_names_gp=grid::gpar(fontsize=6), 
        col=col,
        border=FALSE, 
        cluster_rows=FALSE, 
        cluster_columns=FALSE,
        show_row_names=TRUE,
        show_column_names=TRUE, 
        width=ncol(mat)*unit(2.5, "mm"), 
        height=ncol(mat)*unit(2.5, "mm"), 
        rect_gp=gpar(col="black", lwd=0.5), 
        heatmap_legend_param=list(title_gp=gpar(fontsize=8, fontface="bold"), labels_gp=gpar(fontsize=6))

        ) %>% as.ggplot()
    
    return(hm)
    
    
}
