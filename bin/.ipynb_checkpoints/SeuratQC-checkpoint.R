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

###############
### dplot_1 ###
###############
dplot_theme <- theme(
    
    aspect.ratio=1, legend.position="none",
    axis.title=element_blank(),
    axis.text=element_blank(),
    axis.ticks=element_blank(),
    axis.line=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.border=element_blank()

)

dplot <- function(so, reduction="umap", group_by=NULL, split_by=NULL, label=FALSE, ncol=NULL, legend_position="right", pt_size=0.01, alpha=0.4) {
    
    dplot <- DimPlot(so, reduction=reduction, group.by=group_by, split.by=split_by, label=label, ncol=ncol, pt.size=pt_size, raster=FALSE) & dplot_theme
    if(legend_position!="none") {dplot <- dplot + theme(legend.position=legend_position)}
    dplot[[1]]$layers[[1]]$aes_params$alpha <- alpha
    dplot <- dplot + guides(colour=guide_legend(override.aes=list(alpha=1, size=2), ncol=ncol))
    
    return(dplot)
    
}

###############
### fplot_1 ###
###############
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

fplot <- function(so, reduction="umap", features=NULL, ncol=NULL, legend_position="right", title=NULL, pt_size=0.01, alpha=0.4) {
    
    fplot <- FeaturePlot(so, reduction=reduction, order=TRUE, features=features, pt.size=pt_size) & fplot_theme
    if(legend_position!="none") {fplot <- fplot + theme(legend.position=legend_position)}
    fplot <- fplot + guides(color=guide_colourbar(title=title, title.vjust=5))
    
    return(fplot)
    
}