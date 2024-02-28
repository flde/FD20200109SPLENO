############################
### GSEA from Seurat DEA ###
############################
gsea_from_dea <- function(dea, gene_sets, ident="ident") {
    
    dea$p_val_adj <- ifelse(dea$p_val_adj==0, min(dea$p_val_adj[dea$p_val_adj>0]), dea$p_val_adj)
    dea$sign_log_adj_p_values <- -log10(dea$p_val_adj) * sign(dea$avg_log2FC)
    
    ranks <- dea$sign_log_adj_p_values
    names(ranks) <- dea$gene
    ranks <- ranks[order(ranks)]
    ranks <- rev(ranks)
    
    # Retain only pathways with at least one DEA 
    gene_sets_filter <- lapply(gene_sets, function(x) {sum(dea[dea$p_val_adj <= 0.05, ]$gene %in% x)>=1})
    gene_sets <- gene_sets[unlist(gene_sets_filter)]
    
    # Filter gene set for genes intersect with the DEA gene list 
    gene_sets <- lapply(gene_sets, function(x) {x[x %in% dea$gene]})
    
    # Filter gene sets by number of genes
    gene_sets_filter <- lapply(gene_sets, function(x) {length(x)>=5})
    gene_sets <- gene_sets[unlist(gene_sets_filter)]

    gsea <- fgsea(
        
        pathways=gene_sets,
        stats=ranks,
        nperm=100000, 
        minSize=5,
        maxSize=500
        
    )
    
    gsea$ident <- ident
    
    return(gsea)
    
}

###############################
### GSEA plot for treatment ###
###############################
gsea_plot <- function(gsea, pval_thr, pathway_suffix=NULL, top=20) {
    
    gsea <- as.data.frame(gsea)
    gsea <- na.omit(gsea) 
    
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
    gsea$treatment <- ifelse(sign(gsea$ES)==1, "CpG", "NaCl")
    gsea$treatment <- ifelse(gsea$pval<=pval_thr, gsea$treatment, NA)

    gsea$sign_log_pval_values <- -log10(gsea$pval) * sign(gsea$ES)
    
    # Order hits 
    gsea <- gsea[rev(order(gsea$sign_log_pval_values)), ]
    gsea$pathway <- factor(gsea$pathway, levels=gsea$pathway)

    x_max <- max(abs(gsea$sign_log_pval_values)) + 0.1
    if(x_max<abs(log10(pval_thr))) {x_max <- abs(log10(pval_thr)) + 0.1}

    title <- gsea$ident[1]
    
    plot <- ggplot(gsea, aes(x=abs(log10(pval))*sign(NES), y=pathway, title=title, size=abs(NES), color=treatment)) + 
        geom_point() + 
        ggtitle(title) + xlab("Signed -log10(pval)") + ylab("") +
        xlim(-x_max, x_max) +
        geom_vline(xintercept=abs(log10(pval_thr)), linetype="dashed") +
        geom_vline(xintercept=-1*abs(log10(pval_thr)), linetype="dashed") +         
        scale_color_manual(values=color$treatment, na.value="gray") +
        scale_size(range=c(0, 2)) + 
        guides(
            color=guide_legend(order=1, size=2, keywidth=0.75, keyheight=0.75), 
            size=guide_legend(order=2, title="Abs. (NES)", keywidth=0.75, keyheight=0.75)
        ) + 
        theme(axis.text.y=element_text(size=4, hjust=1, vjust=0.5))
    
    return(plot)
    
}