###############
### library ###
###############

library_load <- suppressMessages(
    
    list(
        
        # Data 
        library(tidyverse), 
        
        # Plotting 
        library(ggplot2), 
        library(patchwork), 
        library(viridis), 
        library(RColorBrewer)
        
    )
    
)

######################
### theme_set_text ###
######################
theme_set_text <- function(size_select=1) {
    
    size <- data.frame(
    
        size_1=c(16, 18, 20),
        size_2=c(10, 12, 16),
        size_3=c(8, 10, 12), 
        size_4=c(6, 6, 8)
        
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

#######################
### density_plot_qc ###
#######################
density_plot_qc <- function(so, title, x, xlab, min, max, fill=time_point, formular=NULL, log10=TRUE, xlim=NULL, nrow=1) {
    
    if(log10) {
        
        plot <- ggplot(so@meta.data, aes(x=log10({{x}}+1), fill={{fill}})) +    
            geom_density(fill="#56B1F7") + 
            geom_vline(aes(xintercept=log10({{min}})), color="red", linetype="longdash") +
            geom_vline(aes(xintercept=log10({{max}})), color="red", linetype="longdash")
        
    } else {
        
        plot <- ggplot(so@meta.data, aes(x={{x}}, fill={{fill}})) +      
            geom_density(fill="#56B1F7") + 
            geom_vline(aes(xintercept={{min}}), color="red", linetype="longdash") +
            geom_vline(aes(xintercept={{max}}), color="red", linetype="longdash")
    }
    
    plot <- plot +
        ggtitle(title) + xlab(xlab) + ylab("Density") +
        facet_wrap(formular, nrow=nrow) + 
        theme(legend.position="bottom", aspect.ratio=1)
    
    if(!is.null(xlim)) {plot <- plot + xlim(xlim)}
    
    return(plot)
    
}

########################
### scattern_plot_qc ###
########################
scattern_plot_qc <- function(so, title, fill, formular=NULL, nrow=1) {
    
    plot <- ggplot(so@meta.data, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA), color={{fill}})) + 
        geom_point(size=3) + 
        ggtitle(title) + ylab("log10(feature count)") + xlab("log10(umi count)") + 
        geom_vline(aes(xintercept=log10(nCount_RNA_min)), color="red", linetype="longdash") +
        geom_vline(aes(xintercept=log10(nCount_RNA_max)), color="red", linetype="longdash") +
        geom_hline(aes(yintercept=log10(nFeature_RNA_min)), color="red", linetype="longdash") + 
        geom_hline(aes(yintercept=log10(nFeature_RNA_max)), color="red", linetype="longdash") + 
        facet_wrap(formular, nrow=nrow) + 
        theme(aspect.ratio=1, legend.position="bottom") + 
        scale_size(guide=guide_legend(direction="vertical"))
    
    return(plot)

}

###################
### box_plot_qc ###
################### 
box_plot_qc <- function(so, y, fill, ylab, ymin, ymax=NULL, formular=NULL) {

    if(is.null(ymax)) {
        ymax=max(so[[deparse(substitute(y))]])
    }
    
    plot <- function(so, title) {
        
        ggplot(so@meta.data, aes(x=sample_name, y={{y}}, color={{fill}})) + 
            geom_jitter(alpha=0.2, shape=16, color="gray") + 
            geom_boxplot(color="#56B1F7", alpha=1.0) + xlab("") + 
            ylim(ymin, ymax) +
            ggtitle(title) + 
            facet_wrap(formular, nrow=1, scales="free_x") + ylab(ylab) +
            theme(
                axis.text.x=element_text(angle=45, vjust=1, hjust=1), 
                legend.position="none"
            )
    }
    
    plot_1 <- plot(so, ylab)
    plot_2 <- plot(subset(so, subset=qc_class=="pass"), paste(ylab, ("(filtered)")))    

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

dplot <- function(so, reduction="umap", group_by=NULL, group_color=NULL, split_by=NULL, label=FALSE, label_size=16, label_box=FALSE, label_color="white", ncol=NULL, legend_position="right", pt_size=2, alpha=1, shape=16, stroke=0, size_select=1, na_color="#7f7f7f", shuffle=TRUE, order=NULL) {
    
    dplot <- DimPlot(so, reduction=reduction, group.by=group_by, split.by=split_by, label=label, label.size=label_size*5/14, label.box=label_box, label.color=label_color, ncol=ncol, pt.size=pt_size, raster=FALSE, na.value=na_color, shuffle=shuffle, order=order) & theme_set_text(size_select=size_select) & dplot_theme 
    if(legend_position!="none") {dplot <- dplot + theme(legend.position=legend_position)}
    dplot[[1]]$layers[[1]]$aes_params$alpha <- alpha
    dplot[[1]]$layers[[1]]$aes_params$shape <- shape


    
    dplot <- dplot + guides(color=guide_legend(override.aes=list(alpha=1, size=4), ncol=ncol, keywidth=0, keyheight=0.75, default.unit="cm"))

    if(!is.null(group_color)) {

        dplot <- dplot + scale_color_manual(values=group_color)
        
    }
    
    
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

fplot <- function(so, reduction="umap", features=NULL, restrict=NULL, order=TRUE, legend_position="right", pt_size=2, alpha=1, shape=16, slot="scale.data", min_cutoff=NA, max_cutoff=NA, min_set=NA, max_set=NA, bar_breaks=waiver(), size_select=1, assay="RNA", option="G") {
    
    DefaultAssay(so)=assay

    scale_select <- 1/size_select

    # Color bar limits 
    if(features %in% rownames(so)) {
        
        color_bar_min <- min(GetAssayData(so, assay=assay, layer=slot)[features, ])
        color_bar_max <- max(GetAssayData(so, assay=assay, layer=slot)[features, ])
        
        suppressMessages(if(slot=="counts") {
            
            color_bar_scale <- scale_color_viridis(option=option, limits=c(0, color_bar_max))
        
        })
        
        suppressMessages(if(slot=="data") {

            color_bar_max <- max(GetAssayData(so, assay=assay, layer=slot)[features, ])
            if(is.na(max_set)) {color_bar_max <- max(GetAssayData(so, assay=assay, layer=slot)[features, ])} else {color_bar_max <- max_set}
            color_bar_max <- ceiling(color_bar_max)
            
            color_bar_scale <- scale_color_viridis(option=option, limits=c(0, color_bar_max), breaks=c(0, color_bar_max/2, color_bar_max), labels=c(0, "", color_bar_max))
        
        }
                        
                        )
        
        suppressMessages(if(slot=="scale.data") {
            
            olor_bar_scale <- scale_color_viridis(option=option, limits=c(color_bar_min, color_bar_max))
            
        })

        so$features <- so@assays[[assay]]@layers[[slot]][rownames(so)==features, ]


        
    } else {
        
        color_bar_min <- ifelse(is.na(min_cutoff), min(na.omit(so[[features]])), min_cutoff)
        color_bar_max <- ifelse(is.na(max_cutoff), max(na.omit(so[[features]])), max_cutoff)
        
        color_bar_scale <- scale_color_viridis(option=option, limits=c(color_bar_min, color_bar_max), breaks=bar_breaks)

        so$features <- so@meta.data[[features]]
        
    }

    # Restrict
    if(!is.null(restrict)) {

        so@meta.data[["features"]][!so@meta.data[[restrict[1]]] %in% restrict[2]] <- NA
    
    }
    
    # Feature plot
    fplot <- FeaturePlot(so, reduction=reduction, order=order, features="features", pt.size=pt_size, min.cutoff=min_cutoff, max.cutoff=max_cutoff, slot=slot) 

    # Always set NA to the back 
    fplot[[1]]$data <- rbind(
    
        fplot[[1]]$data[is.na(fplot[[1]]$data$features), ], 
        fplot[[1]]$data[!is.na(fplot[[1]]$data$features), ] 
        
    )
    
    # Aesthetics 
    fplot <- fplot & theme_set_text(size_select=size_select) & fplot_theme 
    
    if(legend_position!="none") {fplot <- fplot + theme(legend.position=legend_position)}
    
    fplot[[1]]$layers[[1]]$aes_params$alpha <- alpha
    fplot[[1]]$layers[[1]]$aes_params$shape <- shape
    
    suppressMessages(fplot <- fplot + color_bar_scale)
    
    fplot <- fplot + ggtitle(features) + guides(color=guide_colourbar(barwidth=scale_select*1.5, barheight=scale_select*4.0, title="log(CPT)", labels=c(0, 3))) + xlab("") + ylab("")
    
    return(fplot)
    
}