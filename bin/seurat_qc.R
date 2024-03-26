###############
### library ###
###############

library_load <- suppressMessages(
    
    list(
        
        library(patchwork), 
        library(viridis)
        
    )
    
)

######################
### theme_set_text ###
######################
theme_set_text <- function(size_select=1) {
    
    size <- data.frame(
    
        size_1=c(16, 18, 20),
        size_2=c(10, 12, 16),
        size_3=c(8, 10, 12)
        
    )

    size <- size[, size_select]
    
    library(ggplot2, quietly=TRUE)
  
    theme(
        
        panel.spacing=unit(0.2, "lines"),

        plot.title=element_text(size=size[3], face="plain", margin=margin(t=0, r=0, b=10, l=0), color="black"), 
    
        axis.title.y=element_text(size=size[2], face="plain", margin=margin(t=0, r=5, b=0, l=0), angle=90, color="black"), 
        axis.title.x=element_text(size=size[2], face="plain", margin=margin(t=5, r=0, b=0, l=0), color="black"),
        axis.text.y=element_text(size=size[1], face="plain", margin=margin(t=0, r=2, b=0, l=0), color="black"),
        axis.text.x=element_text(size=size[1], face="plain", margin=margin(t=2, r=0, b=0, l=0), color="black"),
      
        axis.ticks=element_line(color="black", size=unit(1/2.141959, "pt")), 
        axis.line=element_line(color="black", size=unit(1/2.141959, "pt")), 
      
        legend.title=element_text(size=size[2], face="plain", color="black"),
        legend.text=element_text(size=size[2], face="plain", color="black"), 
        legend.key=element_rect(fill="transparent", colour="transparent"), 
    
        strip.text=element_text(size=size[2], margin=margin(t=2, r=2, b=2, l=2), face="plain", color="black") 
    
    ) 

}


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
    axis.title.y=element_blank(), 
    axis.title.x=element_blank(),
    axis.text.y=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks=element_blank(),
    axis.line=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.border=element_blank(), 
    legend.key=element_blank()

) 

dplot <- function(so, reduction="umap", group_by=NULL, split_by=NULL, label=FALSE, label_size=16, label_box=FALSE, label_color="white", ncol=NULL, legend_position="right", pt_size=2, alpha=1, shape=16, stroke=0, size_select=1, na_color="grey50") {
    
    dplot <- DimPlot(so, reduction=reduction, group.by=group_by, split.by=split_by, label=label, label.size=label_size*5/14, label.box=label_box, label.color=label_color, ncol=ncol, pt.size=pt_size, raster=FALSE, na.value=na_color) & theme_set_text(size_select=size_select) & dplot_theme 
    if(legend_position!="none") {dplot <- dplot + theme(legend.position=legend_position)}
    dplot[[1]]$layers[[1]]$aes_params$alpha <- alpha
    dplot[[1]]$layers[[1]]$aes_params$shape <- shape
    
    dplot <- dplot + guides(color=guide_legend(override.aes=list(alpha=1, size=4), ncol=ncol, keywidth=0, keyheight=0.75, default.unit="cm"))
    
    
    return(dplot)
    
}

#############
### fplot ###
#############
fplot_theme <- theme(
    
    aspect.ratio=1,
    axis.title.y=element_blank(), 
    axis.title.x=element_blank(),
    axis.text.y=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks=element_blank(),
    axis.line=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.border=element_blank()
    
)

fplot <- function(so, subset=NULL, reduction="umap", features=NULL, legend_position="right", pt_size=2, alpha=1, shape=16, slot="counts", min_cutoff=NA, size_select=1) {
    
    # Color bar limits 
    if(features %in% rownames(so)) {
        
        color_bar_min <- min(GetAssayData(so, assay="RNA", slot=slot)[features, ])
        color_bar_max <- max(GetAssayData(so, assay="RNA", slot=slot)[features, ])
        
        suppressMessages(if(slot=="counts") {color_bar_scale <- scale_color_viridis(option="F", limits=c(0, color_bar_max))})
        suppressMessages(if(slot=="data") {color_bar_scale <- scale_color_viridis(option="F", limits=c(0, color_bar_max))})
        suppressMessages(if(slot=="scale.data") {color_bar_scale <- scale_color_viridis(option="F", limits=c(color_bar_min, color_bar_max))})
        
    }
    
    if(features %in% colnames(so@meta.data)) {
        
        color_bar_min <- min(na.omit(so[[features]]))
        color_bar_max <- max(na.omit(so[[features]]))
        
        color_bar_scale <- scale_color_viridis(option="F", limits=c(color_bar_min, color_bar_max))
        
    }
    
    # Subset 
    if(!is.null(subset)) {so <- so[, so[[subset[1]]]==subset[2]]}
    
    # Feature plot
    fplot <- FeaturePlot(so, reduction=reduction, order=TRUE, features=features, pt.size=pt_size, min.cutoff=min_cutoff, slot=slot) 
    
    # Aesthetics 
    fplot <- fplot & theme_set_text(size_select=size_select) & fplot_theme 
    
    if(legend_position!="none") {fplot <- fplot + theme(legend.position=legend_position)}
    
    fplot[[1]]$layers[[1]]$aes_params$alpha <- alpha
    fplot[[1]]$layers[[1]]$aes_params$shape <- shape
    
    suppressMessages(fplot <- fplot + color_bar_scale)
    
    fplot <- fplot + guides(color=guide_colourbar(barwidth=0.75,barheight=5.5)) + xlab("") + ylab("")
    
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

##################
### dp_feature ###
##################
dp_feature <- function(so, features, split=NULL, group_by, title=NULL, scale=TRUE, assay="RNA", range_min=0, range_max=5) {
    
    # Extract counts 
    mat <- GetAssayData(so, assay=assay, slot="data")[features, ] %>% as.data.frame() %>% add_rownames(var="gene")
    mat <- reshape2::melt(mat, id.vars="gene", value.name="expression", variable.name="cell_id")
    
    # Combine counts with grouping var and compute variables for plotting 
    mat <- dplyr::left_join(mat, so@meta.data[c("cell_id", group_by)], by="cell_id") %>%    
        dplyr::group_by_at(c(group_by, "gene")) %>% 
        dplyr::summarise(mean_expression=mean(expression, na.rm=TRUE), ratio=sum(expression>0)/n(), cell_count=n()) %>%
        dplyr::ungroup() %>% 
        dplyr::mutate(mean_expression=ifelse(is.nan(mean_expression), 0, mean_expression), ratio=ifelse(ratio==0, NA, ratio))
    
    if(scale) {
        
        mat <- dplyr::group_by_at(mat, "gene") %>% dplyr::mutate(mean_expression=scale(mean_expression)) %>% dplyr::ungroup()
        mat <- na.omit(mat)
        
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
        geom_point(aes(size=ratio, fill=mean_expression), alpha=1, shape=21, stroke=0.2) + 
        xlab("") + ylab("") + ggtitle(title) + 
        scale_size_continuous(name="Fraction", limits=c(0, 1), breaks=c(0, 0.25, 0.5, 0.75, 1.0), range=c(range_min, range_max)) + 
        scale_fill_continuous(name="log10(CPT)", low="white", high=rev(brewer.pal(11,"RdBu"))[11]) + 
        theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1)) + 
        theme_global_set()
    
    if(scale) {
        
        limit_scale <- ceiling(max(abs(mat$mean_expression)))
        
        dp <- dp + scale_fill_gradient2(name="z-score", low=rev(brewer.pal(11,"RdBu"))[1], high=rev(brewer.pal(11,"RdBu"))[11], breaks=seq(-limit_scale, limit_scale, 2), limits = c(-limit_scale, limit_scale))
        
    }
    
    if(!is.null(split)) {
        
        dp <- dp + facet_grid(~split, scales="free", space="free")
        
        
    }

    return(dp)
    
}
