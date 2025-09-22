####################
### Library load ###
####################
library_load <- suppressMessages(
    
    suppressWarnings(
        
        list(

            library(msigdbr), 
            library(fgsea)

        )
    )
)

#####################
### GSEA from DEA ###
#####################
gsea_from_dea <- function(dea, gene_set, gene_name=NULL) {
    
    # Set gene names 
    if(!is.null(gene_name)) {dea["gene_name"] <- dea[gene_name]}
    
    # Check format from DEA  
    if(c("logFC") %in% colnames(dea)) {colnames(dea)[colnames(dea)=="logFC"] <- "avg_log2FC"}
    if(c("adj.P.Val") %in% colnames(dea)) {colnames(dea)[colnames(dea)=="adj.P.Val"] <- "p_val_adj"}
    
    # Make ranks 
    dea$p_val_adj <- ifelse(dea$p_val_adj==0, min(dea$p_val_adj[dea$p_val_adj>0]), dea$p_val_adj)
    dea$sign_log_adj_p_values <- -log10(dea$p_val_adj) * sign(dea$avg_log2FC)
    
    ranks <- dea$sign_log_adj_p_values
    names(ranks) <- dea$gene_name
    ranks <- ranks[order(ranks)]
    ranks <- rev(ranks)
    
    # Retain only pathways that overlap with dea lsit
    gene_set_filter <- lapply(gene_set, function(x) {sum(dea$gene %in% x)>=1})
    gene_set <- gene_set[unlist(gene_set_filter)]

    gsea <- fgsea(
        
        pathways=gene_set,
        stats=ranks,
        nperm=100000, 
        minSize=5,
        maxSize=500
        
    )
    
    return(gsea)
    
}

###########################
### GSEA plot for group ###
###########################
gsea_plot <- function(gsea, pval_thr, title=NULL, color=c(RColorBrewer::brewer.pal(8, "Set1")[2], RColorBrewer::brewer.pal(8, "Set1")[1]), color_names=c("Neg", "Pos"), size_range=5, pathway_suffix=NULL, top=20) {
    
    # Set GSEA data frame 
    gsea <- as.data.frame(gsea)
    gsea <- na.omit(gsea) 
    
    # Set color names 
    color <- setNames(color, color_names)
    
    # Fix pathway names
    if(!is.null(pathway_suffix)) {gsea$pathway <- gsub(pathway_suffix, "", gsea$pathway)}
    gsea$pathway <- gsub("_", " ", gsea$pathway)
    
    # Filter hits 
    gsea_up <- gsea[sign(gsea$NES)==+1, ]
    gsea_down <- gsea[sign(gsea$NES)==-1, ]
    
    gsea_up <- gsea_up[order(gsea_up$pval), ][1:top, ]
    gsea_down <- gsea_down[order(gsea_down$pval), ][1:top, ]
    
    gsea <- rbind(gsea_up, gsea_down)
    gsea <- na.omit(gsea)
    gsea <- distinct(gsea)
    
    # Add color 
    gsea$color <- ifelse(sign(gsea$ES)==-1, names(color)[1], names(color)[2])
    gsea$color <- ifelse(gsea$pval<=pval_thr, gsea$color, NA)

    gsea$sign_log_pval_values <- -log10(gsea$pval) * sign(gsea$ES)
    
    # Order hits 
    gsea <- gsea[rev(order(gsea$sign_log_pval_values)), ]
    gsea$pathway <- factor(gsea$pathway, levels=gsea$pathway)

    x_max <- max(abs(gsea$sign_log_pval_values)) + 0.1
    if(x_max<abs(log10(pval_thr))) {x_max <- abs(log10(pval_thr)) + 0.1}
    
    # Plot 
    plot <- ggplot(gsea, aes(x=sign_log_pval_values, y=pathway, color=color)) + 
        
        geom_vline(xintercept=-log10(pval_thr), linetype="dashed") + 
        geom_vline(xintercept=log10(pval_thr), linetype="dashed") +
    
        geom_point(aes(size=abs(NES))) +

        ggtitle(title) +
        xlab("Signed -log10 adj. p-value") + ylab("") + 
        xlim(-x_max, x_max) + 
        scale_color_manual(values=color, na.value="black", drop=FALSE) +
        scale_size(range=c(0, size_range)) + 
        guides(
            
            color=guide_legend(order=1, title="Group", size=2, keywidth=0.75, keyheight=0.75), 
            size=guide_legend(order=2, title="Abs. (NES)", keywidth=0.75, keyheight=0.75)
            
        ) +
    
        theme(
            
            legend.position="right", 
            legend.justification="top", 
            axis.text.y=element_text(hjust=1, vjust=0.5)
            
        )
    
    return(plot)
    
}