


##################
### dp_feature ###
##################
dp_feature=function(so, features, split=NULL, group_by, title=NULL, scale=TRUE, assay="RNA", range_min=0, range_max=6) {
    
    # Extract counts 
    mat <- GetAssayData(so, assay=assay, slot="data")[features, ] %>% as.data.frame() %>% add_rownames(var="gene")
    mat <- reshape2::melt(mat, id.vars="gene", value.name="expression", variable.name="cell_id")
    
    # Combine counts with grouping var and compute variables for plotting 
    mat <- dplyr::left_join(mat, so@meta.data[c("cell_id", group_by)], by="cell_id") %>%
        dplyr::mutate(expression=ifelse(expression==0, NA, expression))  %>%      
        dplyr::group_by_at(c(group_by, "gene")) %>% 
        dplyr::summarise(mean_expression=mean(expression, na.rm=TRUE), ratio=sum(!is.na(expression))/n(), cell_count=n()) %>%
        dplyr::ungroup() %>% 
        dplyr::mutate(mean_expression=ifelse(is.nan(mean_expression), 0, mean_expression), ratio=ifelse(ratio==0, NA, ratio))
    
    if(scale) {
        
        mat <- dplyr::group_by_at(mat, "gene") %>% 
            dplyr::mutate(mean_expression=scale(mean_expression))
        
    }
    
    # Add gene split if available 
    if(!is.null(split)) {
        
        mat <- dplyr::left_join(mat, data.frame(gene=features, split=split), by="gene")
        
    }
    
    # Order genes for plotting 
    mat$gene <- factor(mat$gene, levels=features)
    
    # Order group for plotting 
    if(is.numeric(so@meta.data[[group_by]])) {
        
        mat[[group_by]] <- as.numeric(mat[[group_by]])
        mat[[group_by]] <- factor(mat[[group_by]], levels=sort(unique(mat[[group_by]])))
    
    }
    
    dp <- ggplot(mat, aes_string(x="gene", y=group_by)) + 
        geom_point(aes(size=ratio, fill=mean_expression), alpha=1, shape=21) + 
        xlab("") + ylab("") + ggtitle(title) + 
        scale_size_continuous(name="Fraction", limits=c(0, 1), breaks=c(0, 0.25, 0.5, 0.75, 1.0), range=c(range_min, range_max)) + 
        scale_fill_continuous(name="log10(CPT)", low="white", high=rev(brewer.pal(11,"RdBu"))[11]) + 
        theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1)) + 
        theme_global_set()
    
    if(scale) {
        
        dp <- dp + scale_fill_gradient2(name="z-score", low=rev(brewer.pal(11,"RdBu"))[1], high=rev(brewer.pal(11,"RdBu"))[11])
        
    }
    
    if(!is.null(split)) {
        
        dp <- dp + facet_grid(~split, scales="free", space="free")
        
        
    }

    return(dp)
    
}



######################
### Heatmap of DEA ###
######################
hm_dea <- function(dea, so, p_val=0.05, top=10, pct_2=NULL, title=NULL, column_name, column_order, filter_anno=NULL, width=0.05, column_title_rot=0, conserved=FALSE, re_map=FALSE) {
    
    # Filter by annotation 
    if(!is.null(filter_anno)) {
        
        dea <- lapply(dea, function(x) {x[x[, filter_anno], ]})
    }
    
    if(conserved) {
        
        # Compute averages across groups 
        dea <- lapply(dea, function(x) {x$p_val_adj <- x$minimump_p_val; return(x)})
        dea <- lapply(dea, function(x) {x$avg_log2FC <- rowMeans(x[, c("NaCl_avg_log2FC", "CpG_avg_log2FC")]); return(x)})
        dea <- lapply(dea, function(x) {x$pct.1 <- rowMeans(x[, c("NaCl_pct.1", "CpG_pct.1")]); return(x)})
        dea <- lapply(dea, function(x) {x$pct.2 <- rowMeans(x[, c("NaCl_pct.2", "CpG_pct.2")]); return(x)})
        
        # Filter by avg_log2FC sign 
        dea <- lapply(dea, function(x) {x[sign(x$NaCl_avg_log2FC)==sign(x$CpG_avg_log2FC), ]})
        
    }
    
    # Filter by p-value 
    dea <- lapply(dea, function(x) {x[x$p_val_adj <= p_val, ]})
    
    # Filter by pct.2 
    if(!is.null(pct_2)) {
        
        dea <- lapply(dea, function(x) {x[x$pct.2<=pct_2, ]})
        
    }
    
    # Order by avg_log2FC
    dea <- lapply(dea, function(x) {x[sign(x$avg_log2FC)==1, ]})
    dea <- lapply(dea, function(x) {x[order(-x$avg_log2FC), ]})

    # Filter by top hits 
    dea <- lapply(dea, function(x) {x[1:top, ]})
    
    # Merge dea results 
    dea <- do.call("rbind", dea)
    dea <- na.omit(dea)
    
    # Re-map cell type names 
    if(re_map) {
        
        dea$cell_type <- NULL
        dea <- dplyr::left_join(dea, re_map_df, by="ident")
        
    }
    
    # Gent normalized count matrix and subset by genes from dea
    mat <- GetAssayData(so, assay="RNA", slot="data", features=dea$gene) %>% as.matrix()
    mat <- mat[rownames(mat) %in% dea$gene, ]
    
    # Scale 
    mat <- t(scale(t(mat)))

    # Order rows by genes
    mat <- mat[dea$gene, ]
    
    
    # Select column and order from Seurat meta data and order mat
    so@meta.data[[column_name]] <- factor(so@meta.data[[column_name]], levels=column_order)
    column_order_df <- so@meta.data
    column_order_df <- column_order_df[order(column_order_df[[column_name]]), ]
    column_order_df <- column_order_df[column_name]
    
    mat <- mat[, rownames(column_order_df)]
    
    # Top annotaiton based on column_order
    top_annotation=HeatmapAnnotation(top_annotation=column_order)

    # Breaks 
    breaks <- seq(-quantile(abs(mat), 0.99), quantile(abs(mat), 0.99), length.out=11)
    breaks <- c(breaks[1:4], -0.1, 0, 0.1, breaks[8:11])
    
    color <- c(rev(brewer.pal(11,"RdBu"))[1:4], "#ffffff", "#ffffff", "#ffffff", rev(brewer.pal(11,"RdBu"))[8:11])
    hm <- ComplexHeatmap::Heatmap(
    
        matrix=mat,
        name="z-score", 
        col=colorRamp2(breaks, color), 
        width=ncol(mat)*unit(width, "mm"), 
        height=nrow(mat)*unit(3, "mm"), 
        row_split=factor(dea$cell_type, levels=column_order),
        border=TRUE, 
        cluster_rows=FALSE, 
        cluster_columns=FALSE,
        show_row_names=TRUE,
        show_column_names=FALSE, 
        row_gap=unit(0, "mm"), 
        column_gap=unit(0, "mm"), 
        row_names_gp=gpar(fontsize=8),
        column_title_rot=column_title_rot, 
        row_title_rot=90, 
        column_split=column_order_df[[column_name]]

    ) %>% as.ggplot() + ggtitle(title) + theme(

                plot.title=element_text(size=22, face="bold", margin=margin(t=0, r=0, b=10, l=0), hjust=0.5), 
                panel.grid=element_blank()

            )
    
    return(hm)

}

##############################
### Re-map cell type names ###
##############################
re_map_df <- data.frame(
    
    ident=c(
        
        "17_2", 
        "15_1", 
        "15_0",
        "15_3",
        "17_0",
        "15_2", 
        "17_1", 
        "2_1", 
        "2_2", 
        "2_0", 
        "13", 
        "5", 
        "12", 
        "11", 
        "9", 
        "7", 
        "3", 
        "0", 
        "1", 
        "4", 
        "6", 
        "8",
        "10",
        "14",
        "16", 
        "18", 
        "19", 
        "20"
        
    ), 
    
    cell_type=c(
        
        "Neutrophil (lineage)", 
        "Basophil (lineage)", 
        "BMCP", 
        "MDP", 
        "GMP",
        "CMP", 
        "HSC", 
        "MEP (1)", 
        "MEP (2)", 
        "MEP (3)", 
        "Pro-Erythroblast (1)", 
        "Pro-Erythroblast (2)", 
        "Pro-Erythroblast (3)", 
        "Pro-Erythroblast (4)", 
        "Erythroblast (1)", 
        "Erythroblast (2)", 
        "Erythroblast (3)", 
        "Erythroblast (4)", 
        "Erythroblast (5)", 
        "Mo (1)", 
        "cDC2 (1)", 
        "Mo (2)", 
        "cDC1", 
        "Mo (3)",
        "Mo (4)", 
        "cDC2 (2)", 
        "RPM", 
        "MoRPM"
        
    )
)