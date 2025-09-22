###############
### library ###
###############

library_load <- suppressMessages(
    
    list(
        
        library(patchwork), 
        library(viridis)
        
    )
    
)

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

dplot <- function(so, reduction="umap", group_by=NULL, split_by=NULL, label=FALSE, label_size=16, label_box=FALSE, label_color="black", ncol=NULL, legend_position="right", pt_size=2, alpha=1, shape=16, stroke=0, size_select=1, na_color="grey50", repel=TRUE, label_box_parse=FALSE, label_text_color="black", shuffle=NULL, order=NULL) {
    
    dplot <- DimPlot(so, reduction=reduction, group.by=group_by, split.by=split_by, label=label, label.size=label_size*5/14, label.box=label_box, label.color=label_color, ncol=ncol, pt.size=pt_size, raster=FALSE, na.value=na_color, repel=repel, shuffle=shuffle, order=order) & theme_set_text(size_select=size_select) 
    dplot <- dplot & dplot_theme 
    
    if(legend_position!="none") {dplot <- dplot + theme(legend.position=legend_position)}
    dplot[[1]]$layers[[1]]$aes_params$alpha <- alpha
    dplot[[1]]$layers[[1]]$aes_params$shape <- shape
    
    if(label_box_parse) {dplot$layers[[2]]$geom_params$parse=TRUE}
    if(label_box) {dplot$layers[[2]]$aes_params$colour=label_text_color}
    
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

fplot <- function(so, subset=NULL, reduction="umap", features=NULL, legend_position="right", pt_size=2, alpha=1, shape=16, slot="data", min_cutoff=NA, size_select=1, viridis_option="F") {
    
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
        
        color_bar_scale <- scale_color_viridis(option=viridis_option, limits=c(color_bar_min, color_bar_max))
        
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