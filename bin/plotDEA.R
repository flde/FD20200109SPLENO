##########################################
# Dotplot for FindConservedMarker output #
##########################################
conserved_markers_dp <- function(data, title="Dotplot", p_val_min=0.01, n_top=20, filter=NULL) {
    
    # Load packages 
    require(magrittr)
    require(tidyverse)
    require(ggplot2)
    
    # Grouping names 
    group_1 <- strsplit(colnames(data)[1], "_")[[1]][[1]]
    group_2 <- strsplit(colnames(data)[6], "_")[[1]][[1]]
    
    # Filter data 
    if(filter=="pos_only") {data <- data[sign(data[2])==1 & sign(data[7])==1, ]}

    # Get data  
    data <- dplyr::filter(data, max_pval<=p_val_min) %>% 
        dplyr::arrange(max_pval) %>%
        dplyr::slice(1:n_top) %>% 
        dplyr::mutate(gene=rownames(.))
    


    col_names <- c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "max_pval", "minimump_p_val", "gene", "group")

    data_1 <- cbind(data[1:5], data[11:13]) %>% dplyr::mutate(group=group_1) %>% magrittr::set_colnames(col_names)
    data_2 <- cbind(data[6:10], data[11:13]) %>% dplyr::mutate(group=group_2) %>% magrittr::set_colnames(col_names)

    col_names <- c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "max_pval", "minimump_p_val")

    data <- rbind(data_1, data_2)
    
    # If Inf avg_log2FC to maximum 
    data <- dplyr::mutate(data, avg_log2FC=ifelse(is.infinite(avg_log2FC), max(abs(data$avg_log2FC[!is.infinite(data$avg_log2FC)])), data$avg_log2FC))
    
    # Set group factor level 
    data <- dplyr::mutate(data, group=factor(group, levels=c(group_1, group_2)))
    
    # Plot 
    dp <- ggplot(data, aes(x=group, y=gene)) + 
        geom_point(aes(size=pct.1, color=avg_log2FC), alpha=1)  +
        xlab("") + ylab("") + ggtitle(title) +
        scale_color_continuous(name="Log2FC", limits=c(-max(abs(data$avg_log2FC)), max(abs(data$avg_log2FC))), type="viridis") + 
        scale_size_continuous(name="Fraction", limits=c(0, 1), breaks=c(0, 0.25, 0.5, 0.75, 1.0), range=c(0, 6)) 
    
    return(dp)
 
}

####################################
### Volcano plot for FindMarkers ###
####################################
dea_volcano_plot <- function(deg, log2_thold=1, adjpvalue_thold=0.05, title=NULL) {
    
    #' Volcano plot from Seurat FindMarker differently expressed genes output data frame 
    #' 
    #' @so Seurat FindMarker output 
    #' @log2_thold log2 threshold 
    #' @adjpvalue_thold Adjusted p-value threshold 
    #' @title Plot title 
    
    # Set rownames to genes
    if("gene" %in% colnames(deg)) {rownames(deg) <- deg$gene}
    
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
    
    # Set labels 
    deg$label <- ifelse(rownames(deg) %in% c(pos_labels, neg_labels), rownames(deg), NA)
    
    # Plot
    library(ggrepel, quietly=TRUE)
    volcano_plot <- ggplot(deg, aes(x=avg_log2FC, y=-log10(p_val_adj), colour=deg$color, label=label)) + 
        geom_point() + 
        geom_vline(aes(xintercept=log2_thold, color="black", linetype="longdash")) +
        geom_vline(aes(xintercept=-log2_thold, color="black", linetype="longdash")) +
        geom_hline(aes(yintercept=-log10(adjpvalue_thold), color="black", linetype="longdash")) + 
        geom_label_repel(segment.color="grey50", force=10, force_pull=1, max.overlaps=getOption("ggrepel.max.overlaps", default=100)) + 
        xlim(-max(abs(deg$avg_log2FC)), max(abs(deg$avg_log2FC))) +  
        ggtitle(title) + xlab("average log2FC") + ylab("-log10(adj. p-value)") + 
        scale_colour_manual(values=color) + 

        theme(aspect.ratio=1, legend.position="none")
    
    return(volcano_plot)
}