library_load <- suppressMessages(
    
    list(
        
        # Data 
        library(tidyverse), 
        
        # Plotting 
        library(ggplot2), 
        library(patchwork), 
        library(RColorBrewer)
        
    )
    
)

##################
### dp_feature ###
##################
dp_feature <- function(so, features, split=NULL, split_order=NULL, group_by, group_by_order=NULL, title=NULL, scale=TRUE, assay="RNA", range_min=0, range_max=5) {
    
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
        
        mat <- dplyr::group_by(mat, gene) %>% 
            # dplyr::mutate(mean_expression=(mean_expression-min(mean_expression))) %>% 
            dplyr::mutate(mean_expression=mean_expression/max(mean_expression)) %>% 
            dplyr::ungroup() %>% na.omit()
        
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
    
    if(!is.null(group_by_order)) {
        
        mat[[group_by]] <- factor(mat[[group_by]], levels=group_by_order)
        
    }
    
    # Order split for plotting 
    if(!is.null(split_order)) {
        
        mat[["split"]] <- factor(mat[["split"]], levels=split_order)
        
    }
    
    dp <- ggplot(mat, aes_string(x="gene", y=group_by)) + 
        geom_point(aes(size=ratio, fill=mean_expression), alpha=1, shape=21, stroke=0.2) + 
        xlab("") + ylab("") + ggtitle(title) + 
        scale_size_continuous(name="Fraction", limits=c(0, 1), breaks=c(0, 0.25, 0.5, 0.75, 1.0), range=c(range_min, range_max)) + 
        scale_fill_continuous(name="log(CPM)", low="white", high=rev(brewer.pal(11,"RdBu"))[11])
    
    # Set scale legend 
    if(scale) {
        
        dp <- dp + scale_fill_gradient2(name="Scaled", low=rev(brewer.pal(11,"RdBu"))[1], high=rev(brewer.pal(11,"RdBu"))[11], breaks=c(0.0, 0.5, 1.0), labels=c(0.0, 0.5, 1.0))
        
    }
    
    # Set split
    if(!is.null(split)) {
        
        dp <- dp + facet_grid(~split, scales="free", space="free")
        
    }
    
    # Add theme 
    dp <- dp + 

        theme(

            axis.text.x=element_text(angle=45, vjust=1, hjust=1, face="italic"),
            axis.text.y=element_text(angle=0, vjust=0.5, hjust=1), 
            strip.text.x=element_text(angle=90, vjust=0.5, hjust=0), 
            strip.background=element_rect(fill="transparent", color=NA), 
            legend.key.size=unit(0.25, 'cm'), 
            legend.key.height=unit(0.25, 'cm'), 
            legend.key.width=unit(0.25, 'cm')
            
        )

    return(dp)
    
}