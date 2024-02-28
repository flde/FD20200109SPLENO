#########################
### Heat map net slot ###
#########################
hm_net <- function(cellchat, slot, title, color_max=RColorBrewer::brewer.pal(9, "OrRd")[9]) {
    
    mat <- cellchat@net[[slot]]
    
    if (slot=="count") {name="Count"}
    if (slot=="weight") {name="Weight"}
    
    # Color and Breaks 
    if(name=="Count") {
        
        breaks <- seq(0, 10*ceiling(max(mat)/10), 1)
        color <- colorRampPalette(c("#ffffff", color_max))(length(breaks))
        
    } else if(name=="Weight") {
        
        breaks <- seq(0, 0.1*ceiling(max(mat)/0.1), 0.1)
        color <- colorRampPalette(c("#ffffff", color_max))(length(breaks))
        
    } 

    hm <- ComplexHeatmap::Heatmap(
    
        matrix=mat,
        name=name, 
        col=colorRamp2(breaks, color), 
        column_title_gp=gpar(fontsize=10, fontface="bold"), 
        column_names_gp =grid::gpar(fontsize=6), 
        row_names_gp=grid::gpar(fontsize=6),
        border=TRUE, 
        cluster_rows=FALSE, 
        cluster_columns=FALSE,
        show_row_names=TRUE,
        show_column_names=TRUE, 
        row_gap=unit(1, "mm"), 
        column_gap=unit(1, "mm"), 
        column_title=title, 
        row_title_rot=90, 
        width=ncol(mat)*unit(2.5, "mm"), 
        height=ncol(mat)*unit(2.5, "mm"), 
        rect_gp=gpar(col="black", lwd=0.5), 
        heatmap_legend_param=list(title_gp=gpar(fontsize=8, fontface="bold"), labels_gp=gpar(fontsize=6))

    ) %>% as.ggplot() + theme(

                plot.title=element_text(size=22, face="bold", margin=margin(t=0, r=0, b=0, l=0), hjust=0.5), 
                panel.grid=element_blank()

            )

    return(hm)
}