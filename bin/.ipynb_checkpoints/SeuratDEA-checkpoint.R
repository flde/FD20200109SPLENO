library_load <- suppressMessages(
    
    list(
        
        # Seurat 
        library(Seurat), 
        
        # Plotting 
        library(ComplexHeatmap)
    )
)

library_load <- suppressMessages(
    
    list(
        
        # Plotting 
        library(ComplexHeatmap), 
        library(RColorBrewer)
    )
)

#############################
### average_expression_hm ###
#############################
average_expression_hm <- function(deg, av_exp, top=25, avg_log2FC="pos", cluster_rows=TRUE, scale="row") {
    
    # Set gene column to rownames if present 
    if("gene" %in% colnames(deg)) {rownames(deg) <- deg$gene}
    
    # Subset differential gene expression matrix by top hits 
    if(avg_log2FC=="pos") {deg <- deg[deg$avg_log2FC > 0, ]}
    if(avg_log2FC=="neg") {deg <- deg[deg$avg_log2FC < 0, ]}
    
    deg <- deg[1:top, ]

    # Subset average expression matrix by top hit features from deg
    av_exp <- av_exp[rownames(av_exp) %in% rownames(deg), ]
    av_exp <- av_exp[rownames(deg), ]
    
    color_ramp <- viridis(100, option="inferno")
    if(scale=="none") {name<-"Exp."} else {name<-"z-score"}
    
    hm <- grid.grabExpr(
        draw(
            ComplexHeatmap::pheatmap(
                mat=av_exp,
                main=paste("Cluster", deg$cluster[1]),
                name=name, 
                fontsize_row=10,
                scale=scale,
                cluster_rows=cluster_rows,
                cluster_cols=FALSE,
                treeheight_row=0, 
                cellwidth=10, 
                cellheight=10, 
                show_rownames=TRUE,
                show_colnames=TRUE,
                color=color_ramp,
                border_color=NA
            )
        )
    )
    
    return(hm)
    
}
