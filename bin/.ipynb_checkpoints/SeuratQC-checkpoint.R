library_load <- suppressMessages(
    
    list(
        # Seurat 
        library(patchwork)
    )
)

####################
### rank_plot_qc ###
####################
rank_plot_qc <- function(so, color_color, formular) {
    
    plot <- ggplot(so@meta.data, aes(x=log10(nCount_RNA_rank), y=log10(nCount_RNA), color=cellranger_class)) + 
        geom_point() + 
        scale_color_manual(values=color_color) +
        ggtitle("Barcode rank plot") +
        xlab("log10(cell barcode rank)") + ylab("log10(cell UMI counts)") + 
        facet_grid(formular) + 
        theme(aspect.ratio=1, legend.position="bottom")
    
    return(plot)
    
}

#######################
### density_plot_qc ###
#######################
density_plot_qc <- function(so, title, x, xlab, min, max, fill=tissue, fill_color, formular=tissue~treatment, log10=TRUE) {
    
    if(log10) {
        
        plot <- ggplot(so@meta.data, aes(x=log10({{x}}+1), fill={{fill}})) + 
            geom_density() +
            geom_vline(aes(xintercept=log10({{min}})), color="red", linetype="longdash") +
            geom_vline(aes(xintercept=log10({{max}})), color="red", linetype="longdash")
        
    } else {
        
        plot <- ggplot(so@meta.data, aes(x={{x}}, fill={{fill}})) +  
            geom_density() +
            geom_vline(aes(xintercept={{min}}), color="red", linetype="longdash") +
            geom_vline(aes(xintercept={{max}}), color="red", linetype="longdash")
    }
    
    plot <- plot +
        ggtitle(title) + xlab(xlab) + ylab("Density") +
        scale_fill_manual(values=fill_color) +
        facet_grid(formular) + 
        theme(legend.position="bottom", aspect.ratio=1)
    
    return(plot)
    
}

########################
### scattern_plot_qc ###
########################
scattern_plot_qc <- function(so, title, color) {
    
    plot <- ggplot(so@meta.data, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA), color={{color}})) + 
        geom_point() + ggtitle(title) + ylab("log10(feature count)") + xlab("log10(umi count)") + 
        geom_vline(aes(xintercept=log10(nCount_RNA_min)), color="red", linetype="longdash") +
        geom_vline(aes(xintercept=log10(nCount_RNA_max)), color="red", linetype="longdash") +
        geom_hline(aes(yintercept=log10(nFeature_RNA_min)), color="red", linetype="longdash") + 
        geom_hline(aes(yintercept=log10(nFeature_RNA_max)), color="red", linetype="longdash") + 
        facet_grid(tissue~treatment) + theme(aspect.ratio=1, legend.position="bottom") + 
        scale_size(guide=guide_legend(direction="vertical"))
    
    return(plot)

}

###################
### box_plot_qc ###
################### 
box_plot_qc <- function(so, y, fill, ylab, ymin, ymax=NULL) {

    if(is.null(ymax)) {
        
        ymax=max(so[[deparse(substitute(y))]])
    }
    
    plot <- function(so, title) {
        
        ggplot(so@meta.data, aes(x=treatment, y={{y}}, color={{fill}})) + 
            geom_jitter(alpha=0.2, shape=16, color="gray") + 
            geom_boxplot(alpha=1.0) + xlab("") + 
            ylim(ymin, ymax) +
            scale_color_manual(values=color$tissue) + 
            ggtitle(title) + 
            facet_wrap(~tissue, scales="free_x") + ylab(ylab) +
            theme(
                plot.title=element_text(size=12, face="bold", margin=margin(t=0, r=0, b=5, l=0)), 
                axis.text.x=element_text(angle=45, vjust=1, hjust=1), 
                strip.text=element_blank()
            )
    }
    
    plot_1 <- plot(subset(so, subset=cellranger_class == "Cell"), "GEM")
    plot_2 <- plot(subset(so, subset=cellranger_class == "Cell" & qc_class =="pass"), "GEM and QC")    

  return(list(plot_1, plot_2))
    
}

########################
### deg_volcano_plot ###
########################
deg_volcano_plot <- function(deg, log2_thold = 1, adjpvalue_thold = 0.05, title = "Volcano plot") {
    
    #' Volcano plot from Seurat FindMarker differently expressed genes output data frame 
    #' 
    #' @so Seurat FindMarker output 
    #' @log2_thold log2 threshold 
    #' @adjpvalue_thold Adjusted p-value threshold 
    #' @title Plot title 
    
    # Annotate entries significance by log2_thold and adjpvalue_thold
    deg$p_val_adj <- ifelse(deg$p_val_adj == 0, .Machine$double.xmin, deg$p_val_adj)
    deg$sig <- ifelse(abs(deg$avg_log2FC) >= log2_thold & -log10(deg$p_val_adj) >= -log10(adjpvalue_thold), "s", "ns")
    
    # Set color based on significance and direction of deg e.g. positive and negative 
    deg$color <- ifelse(deg$sig == "s" & deg$avg_log2FC > 0, "s_pos", "ns")
    deg$color <- ifelse(deg$sig == "s" & deg$avg_log2FC < 0, "s_neg", deg$color)
    
    color <- c(RColorBrewer::brewer.pal(8, "Set1")[1], "gray", "black", RColorBrewer::brewer.pal(8, "Set1")[2])
    names(color) <- c("s_pos", "ns", "black", "s_neg")
    
    # Create labels based log2FC and p_val_adj
    deg_pos <- deg[deg$avg_log2FC > 0 & deg$sig == "s", ]
    deg_neg <- deg[deg$avg_log2FC < 0 & deg$sig == "s", ]

    pos_labels_log2FC <- deg_pos[rev(order(deg_pos$avg_log2FC)), ][1:10, ] %>% rownames()
    neg_labels_log2FC <- deg_neg[order(deg_neg$avg_log2FC), ][1:10, ] %>% rownames()
    
    pos_labels_p_val_adj <- deg_pos[order(deg_pos$p_val_adj), ][1:10, ] %>% rownames()
    neg_labels_p_val_adj <- deg_neg[order(deg_neg$p_val_adj), ][1:10, ] %>% rownames()
    
    pos_labels <- c(pos_labels_log2FC, pos_labels_p_val_adj)
    neg_labels <- c(neg_labels_log2FC, neg_labels_p_val_adj)

    deg$label <- ifelse(rownames(deg) %in% c(pos_labels, neg_labels), rownames(deg), NA)
    
    # Plot
    library(ggrepel, quietly = TRUE)
    volcano_plot <- ggplot(deg, aes(x = avg_log2FC, y = -log10(p_val_adj), colour = deg$color, label = label)) + 
        geom_point() + 
        geom_vline(aes(xintercept = log2_thold, color = "black", linetype = "longdash")) +
        geom_vline(aes(xintercept = -log2_thold, color = "black", linetype = "longdash")) +
        geom_hline(aes(yintercept = -log10(adjpvalue_thold), color = "black", linetype = "longdash")) + 
        geom_label_repel(segment.color = "grey50", force = 10, force_pull = 1, max.overlaps = getOption("ggrepel.max.overlaps", default = 100)) + 
        xlim(-max(abs(deg$avg_log2FC)), max(abs(deg$avg_log2FC))) +  ylim(-20, 350) +
        ggtitle(title) + xlab("average log2FC") + ylab("-log10(adj. p-value)") + 
        scale_colour_manual(values = color) + 

        theme(aspect.ratio = 1, legend.position = "none")
    
    return(volcano_plot)
}

#####################
### pc_annotation ###
#####################
pc_annotation <- function(so, dims = 25, loadings = 100, ident = ident) {
    
    #' Bar plot of top PC with feature annotation
    #' 
    #' @so Seurat object with PCA slot 
    #' @dims Number Of PC dims
    #' @loadings Number of top gene loadings by abs
    #' @ident ident from so to become ident column in output for plotting tilte. 
    
    
    # Get cell annotation vectors 
    library(msigdbr, quietly = TRUE)
    go_bp <- msigdbr(species = "mouse", category = "C5", subcategory = "BP")

    cc_genes <- go_bp[go_bp$gs_exact_source == "GO:0022402", ]$gene_symbol
    mt_genes <- rownames(so)[grep("^mt-", rownames(so))]
    hb_genes <- rownames(so)[grep("Hba-|Hbb-|Hbq1b|Hbq1a", rownames(so))]
    rb_genes <- rownames(so)[grep("^Rpl|^Rps", rownames(so))]
    
    # Extract PC loadings 
    pca_loadings <- Loadings(so, reduction = "pca")[, 1:dims]
    
    # Annotate PC
    pca_loadings_genes <- list()
    for(pc in colnames(pca_loadings)) {

        df <- data.frame(
            
            pca_loadings_genes = names(rev(sort(abs(pca_loadings[, pc])))[1:loadings]), 
            pc = pc, 
            ident = so[[ident]][1, 1]
            
        )

        df$annotation <- ifelse(df$pca_loadings_genes %in% cc_genes, "Cell cycle", "Other")
        df$annotation <- ifelse(df$pca_loadings_genes %in% mt_genes, "Mitochondrial", df$annotation)
        df$annotation <- ifelse(df$pca_loadings_genes %in% hb_genes, "Hemoglobin", df$annotation)
        df$annotation <- ifelse(df$pca_loadings_genes %in% rb_genes, "Ribosomal", df$annotation)
        
        pca_loadings_genes[[pc]] <- df
    }
    
    # Combine results 
    pca_loadings_genes <- do.call("rbind", pca_loadings_genes)
    
    # Set factor levels 
    pca_loadings_genes$pc <- factor(pca_loadings_genes$pc, levels = paste0("PC_", 1:dims))
    pca_loadings_genes$annotation <- factor(pca_loadings_genes$annotation, levels = rev(c("Mitochondrial", "Hemoglobin", "Ribosomal", "Cell cycle", "Other")))
    
    
    return(pca_loadings_genes)
    
}

##########################
### pc_annotation_plot ###
##########################
pc_annotation_plot <- function(df) {
    
    p <- ggplot(df, aes(x = pc, fill = annotation)) + 
        geom_bar(stat = "count", position = "fill") + 
        scale_fill_manual(values = color$module_class) + 
        ggtitle(df$ident[1]) + xlab("") + ylab("Frequency") + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    return(p)
}

###############
### dplot_1 ###
###############

dplot_1 <- function(so, reduction="umap", cluster) {
    
    dplot_1 <- DimPlot(so, reduction=reduction, group.by=cluster, label=TRUE) & 
        theme(aspect.ratio=1, legend.position="none")
    
    dplot_2 <- DimPlot(so, reduction=reduction, group.by="tissue", label=FALSE) & 
        theme(aspect.ratio=1, legend.position="bottom") & 
        scale_color_manual(values=color$tissue, na.value="dark gray") & 
        guides(color=guide_legend(ncol=3, override.aes=list(size=2)))
    
    dplot_3 <- DimPlot(so, reduction=reduction, group.by="cc_phase_class", label=FALSE) & 
        theme(aspect.ratio=1, legend.position="bottom") & 
        scale_color_manual(values=color$cc_phase_class, na.value="dark gray") & 
        guides(color=guide_legend(ncol=3, override.aes=list(size=2)))
    
    dplot_4 <- DimPlot(so, reduction=reduction, group.by="treatment", label=FALSE) & 
        theme(aspect.ratio=1, legend.position="bottom") & 
        scale_color_manual(values=color$treatment, na.value="dark gray") &
        guides(color=guide_legend(ncol=3, override.aes=list(size=2)))
    
    dplot_5 <- DimPlot(so, reduction=reduction, group.by="main_labels", label=FALSE) & 
        theme(aspect.ratio=1, legend.position="bottom") & 
        scale_color_manual(values=color$main_labels, na.value="dark gray") & 
        guides(color=guide_legend(ncol=3, override.aes=list(size=2)))
    
    dplot_6 <- DimPlot(so, reduction=reduction, group.by="fine_labels", label=FALSE) & 
        theme(aspect.ratio=1, legend.position="bottom") & 
        scale_color_manual(values=color$fine_labels, na.value="dark gray") & 
        guides(color=guide_legend(ncol=3, override.aes=list(size=2)))
    
    # Combine plots 
    dplot <- dplot_1 + dplot_2 + dplot_3 + dplot_4 + dplot_5 + dplot_6 + plot_layout(ncol=3) 
    return(dplot)
    
}

###############
### fplot_1 ###
###############

fplot_1 <- function(so, reduction="umap") {
    
    fplot_1 <- FeaturePlot(so, reduction = reduction, features = "nCount_RNA") & theme(aspect.ratio = 1) & ggtitle("UMI count")
    
    fplot_2 <- FeaturePlot(so, reduction = reduction, features = "nFeature_RNA") & theme(aspect.ratio = 1) & ggtitle("Feature count")
    
    fplot_3 <- FeaturePlot(so, reduction = reduction, features = "pMt_RNA") & theme(aspect.ratio = 1) & ggtitle("Mt [%]")
    
    fplot_4 <- FeaturePlot(so, reduction = reduction, features = "pHb_RNA") & theme(aspect.ratio = 1) & ggtitle("Hb [%]")
    
    fplot_5 <- FeaturePlot(so, reduction = reduction, features = "pRp_RNA") & theme(aspect.ratio = 1) & ggtitle("Rp [%]")
    
    fplot_6 <- FeaturePlot(so, reduction = reduction, features = "msMHCI_RNA1") & theme(aspect.ratio = 1) & ggtitle("MHC I")
    
    fplot_7 <- FeaturePlot(so, reduction = reduction, features = "msMHCII_RNA1") & theme(aspect.ratio = 1) & ggtitle("MHC II")
    
    fplot_8 <- FeaturePlot(so, reduction = reduction, features = "msTcr_RNA1") & theme(aspect.ratio = 1) & ggtitle("TCR")
    
    fplot_9 <- FeaturePlot(so, reduction = reduction, features = "msCd4_RNA1") & theme(aspect.ratio = 1) & ggtitle("CD4")
    
    fplot_10 <- FeaturePlot(so, reduction = reduction, features = "msCd8_RNA1") & theme(aspect.ratio = 1) & ggtitle("CD8")
    
    fplot_11 <- FeaturePlot(so, reduction = reduction, features = "msIgkc_RNA1") & theme(aspect.ratio = 1) & ggtitle("Igkc")
    
    fplot_12 <- FeaturePlot(so, reduction = reduction, features = "msIglc_RNA1") & theme(aspect.ratio = 1) & ggtitle("Iglc")
    
    # Combine plots 
    dplot <- fplot_1 + fplot_2 + fplot_3 + fplot_4 + fplot_5 + fplot_6 + fplot_7 + fplot_8 + fplot_9 + fplot_10 + fplot_11 + fplot_12 + plot_layout(ncol=4)
    
    return(dplot)
}