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
density_plot_qc <- function(so, title, x, xlab, min, max, fill=tissue, fill_color, formular=tissue~treatment+sample_rep, log10=TRUE) {
    
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
        facet_grid(tissue~treatment+sample_rep) + theme(aspect.ratio=1, legend.position="bottom") + 
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
        
        ggplot(so@meta.data, aes(x=sample_rep, y={{y}}, color={{fill}})) + 
            geom_jitter(alpha=0.2, shape=16, color="gray") + 
            geom_boxplot(alpha=1.0) + xlab("") + 
            ylim(ymin, ymax) +
            scale_color_manual(values=color$tissue) + 
            ggtitle(title) + 
            facet_grid(~tissue+treatment, scales="free_x") + ylab(ylab) +
            theme(
                plot.title=element_text(size=12, face="bold", margin=margin(t=0, r=0, b=5, l=0)), 
                axis.text.x=element_text(angle=45, vjust=1, hjust=1)
            )
    }
    
    plot_1 <- plot(subset(so, subset=cellranger_class == "Cell"), "GEM")
    plot_2 <- plot(subset(so, subset=cellranger_class == "Cell" & qc_class =="pass"), "GEM and QC")    

  return(list(plot_1, plot_2))
    
}

#############
### dplot ###
#############
dplot_theme <- theme(
    
    aspect.ratio=1, 
    legend.position="none",
    axis.title=element_blank(),
    axis.text=element_blank(),
    axis.ticks=element_blank(),
    axis.line=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.border=element_blank(), 

    legend.key=element_blank()

)

dplot <- function(so, reduction="umap", group_by=NULL, split_by=NULL, label=FALSE, ncol=NULL, legend_position="right", pt_size=1, alpha=1, shape=16, stroke=0) {
    
    dplot <- DimPlot(so, reduction=reduction, group.by=group_by, split.by=split_by, label=label, ncol=ncol, pt.size=pt_size, raster=FALSE) & dplot_theme
    if(legend_position!="none") {dplot <- dplot + theme(legend.position=legend_position)}
    dplot[[1]]$layers[[1]]$aes_params$alpha <- alpha
    dplot[[1]]$layers[[1]]$aes_params$shape <- shape
    
    dplot <- dplot + guides(color=guide_legend(override.aes=list(alpha=1, size=3), ncol=ncol, keywidth=NULL, keyheight=NULL))
    
    return(dplot)
    
}

#############
### fplot ###
#############
fplot_theme <- theme(
    aspect.ratio=1,
    axis.title=element_blank(),
    axis.text=element_blank(),
    axis.ticks=element_blank(),
    axis.line=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.border=element_blank()
)

fplot <- function(so, subset=NULL, reduction="umap", features=NULL, legend_position="right", pt_size=1, alpha=1, shape=16, slot="scale.data", min_cutoff=NA) {
    
    # Color bar limits 
    if(features %in% rownames(so)) {
        
        color_bar_min <- min(GetAssayData(so, assay="RNA", slot=slot)[features, ])
        color_bar_max <- max(GetAssayData(so, assay="RNA", slot=slot)[features, ])
        
        suppressMessages(if(slot=="counts") {color_bar_scale <- scale_color_viridis(option="F", limits = c(0, color_bar_max))})
        suppressMessages(if(slot=="data") {color_bar_scale <- scale_color_viridis(option="F", limits = c(0, color_bar_max))})
        suppressMessages(if(slot=="scale.data") {color_bar_scale <- scale_color_viridis(option="F", limits = c(color_bar_min, color_bar_max))})
        
    }
    
    if(features %in% colnames(so@meta.data)) {
        
        color_bar_min <- min(na.omit(so[[features]]))
        color_bar_max <- max(na.omit(so[[features]]))
        
        color_bar_scale <- scale_color_viridis(option="F", limits = c(color_bar_min, color_bar_max))
        
    }
    
    # Subset 
    if(!is.null(subset)) {so <- so[, so[[subset[1]]]==subset[2]]}
    
    # Feature plot
    fplot <- FeaturePlot(so, reduction=reduction, order=TRUE, features=features, pt.size=pt_size, min.cutoff=min_cutoff, slot=slot) 
    
    # Aesthetics 
    fplot <- fplot & fplot_theme
    
    if(legend_position!="none") {fplot <- fplot + theme(legend.position=legend_position)}
    
    fplot[[1]]$layers[[1]]$aes_params$alpha <- alpha
    fplot[[1]]$layers[[1]]$aes_params$shape <- shape
    
    suppressMessages(fplot <- fplot + color_bar_scale)
    
    fplot <- fplot + guides(color=guide_colourbar(barwidth=0.75,barheight=5.5))
    
    return(fplot)
    
}

##############
### vjplot ###
##############
vjplot <- function(so, group_by, split_by, genes) {
    
    so <- NormalizeData(so, normalization.method="LogNormalize", scale.factor=10000)
    mat <- GetAssayData(so, assay="RNA", slot="data")
    mat <- as.matrix(mat)
    mat <- as.data.frame(t(mat))
    mat <- mat[, genes, drop=FALSE]
    
    meta <- so@meta.data[, c({group_by}, {split_by}), drop=FALSE]
    
    data <- dplyr::left_join(dplyr::mutate(meta, cell_id=rownames(meta)), dplyr::mutate(mat, cell_id=rownames(mat)), by="cell_id")
    data <- dplyr::select(data, -cell_id)
    
    data <- reshape2::melt(data, id.vars=c({group_by}, {split_by}))
    data <- split(data, f=data$variable)
    
    plot <- function(data) {
        
        y_max <- max(data$value) + 0.25
        
        ggplot(data, aes_string(x=group_by, y="value", fill=split_by, color=split_by)) + 
            geom_violin(alpha=0.5, position=position_dodge(width=0.75), width=1.0, color=NA) +
#             geom_boxplot(aes(outlier.colour=treatment), position=position_dodge(width=1), width=0.1, size=0.25, outlier.size=1.0, outlier.shape=21, outlier.colour=NA) + 
            geom_jitter(position=position_jitterdodge(jitter.width=0.5, dodge.width=0.75), size=1) + 
            xlab("") + ylab("log(CPT+1)") + 
            ylim(0, y_max) + 
            scale_fill_manual(values=color$treatment) + 
            scale_color_manual(values=color$treatment) +
            ggtitle(data$variable) + theme(legend.position="none", aspect.ratio=1) + theme_global_set()
    
    }

    plot_list <- lapply(data, plot)
    
    return(plot_list)
    
}