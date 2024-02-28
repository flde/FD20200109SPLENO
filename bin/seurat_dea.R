####################
### filter_genes ###
####################
feature_select <- function(so, cnt_min=1, cell_min=5, ident_i=NULL, gene_filter=NULL) {
    
    # Subset the count matrix if ident is present
    if(!is.null(ident)) {
        
        so_tmp <- subset(so, idents=ident_i)
        
    } else {
        
        so_tmp <- so
    }
    
    cnt <- GetAssayData(so_tmp, assay="RNA", slot="counts")
    cnt <- cnt[rowSums(cnt>=cnt_min)>=cell_min, ]
    
    if(!is.null(gene_filter)) {
        
        message(paste0("Filter genes: Found ", sum(rownames(cnt) %in% gene_filter), " of ", length(gene_filter)))
        cnt <- cnt[rownames(cnt) %in% gene_filter, ]
        
    }
    
    
    so_tmp <- CreateSeuratObject(counts=cnt, meta.data=so_tmp@meta.data)
    so_tmp <- NormalizeData(so_tmp)
    
    
    
    return(so_tmp)
    
}

###############
### dea_FUN ###
###############
dea_seurat <- function(so, ident="seurat_clusters",  map=NULL, file=NULL, only_pos=TRUE, logfc_threshold=0, min_pct=0, conserved=FALSE, treatment=FALSE, grouping_var=NULL, compute=TRUE, test_use="wilcox", cnt_min=3, cell_min=3, gene_filter=NULL) {
    
    # Load results only if compute FALSE
    if(!compute) {dea <- readRDS(paste0(file, ".rds")); return(dea)}
    
    # Set Ident for comparison 
    so <- SetIdent(so, value=ident)
    
    # Build dummy mapping of idents to cell_types 
    if(is.null(map)) {
        
        map <- data.frame(ident=unique(so@meta.data[[ident]]), cell_type=unique(so@meta.data[[ident]]))
        map[[ident]] <- unique(so@meta.data[[ident]])
        
    }
    
    # DEA
    dea <- list()
    for(n in 1:nrow(map)) {
        
        # Select the column of map to split the Seurat object for DEA by ident 
        i <- map[, ident][n] 
        
        if (conserved) {
            
            message("Run DEA mode: Conserved")
            so_tmp <- feature_select(so, ident=NULL, cnt_min=cnt_min, cell_min=cell_min, gene_filter=gene_filter)
            so_tmp <- SetIdent(so_tmp, value=ident)
            dea_result <- FindConservedMarkers(so_tmp, ident.1=i, logfc.threshold=logfc_threshold, min.pct=min_pct, grouping.var=grouping_var, only.pos=only_pos, subset.ident=NULL, test.use=test_use)
        
        # Treatment DEA
        } else if (treatment) {
            
            message("Rund DEA mode: Treatment")
            so_tmp <- feature_select(so, ident=i, cnt_min=cnt_min, cell_min=cell_min, gene_filter=gene_filter)
            so_tmp <- SetIdent(so_tmp, value=grouping_var)
            dea_result <- FindMarkers(so_tmp, ident.1="CpG", ident.2="NaCl", logfc.threshold=logfc_threshold, min.pct=min_pct, group.by=grouping_var, only.pos=only_pos, subset.ident=NULL, test.use=test_use)
        
        # Marker DEA
        } else {
            
            message("Run DEA mode: Marker")
            so_tmp <- feature_select(so, ident=NULL, cnt_min=cnt_min, cell_min=cell_min, gene_filter=gene_filter)
            so_tmp <- SetIdent(so_tmp, value=ident)
            dea_result <- FindMarkers(so_tmp, ident.1=i, logfc.threshold=logfc_threshold, min.pct=min_pct, group.by=grouping_var, only.pos=only_pos, subset.ident=NULL, test.use=test_use)
            
        }
        
        # Annotate results 
        if(nrow(dea_result)>0) {
            
            # Add meta data to result DEA
            dea_result$gene <- rownames(dea_result)
            dea_result$ident <- i

            # Annotate transcription factors 
            tf <- read.delim("data/annotation/animaltfdb3/Mus_musculus_TF.txt")
            dea_result$transcription_factor <- ifelse(dea_result$gene %in% tf$Symbol, TRUE, FALSE)

            # Annotate receptor and ligands 
            CellChatDB <- CellChat::CellChatDB.mouse
            dea_result$ligand  <- ifelse(dea_result$gene %in% CellChatDB$interaction[["ligand"]], TRUE, FALSE)
            dea_result$receptor  <- ifelse(dea_result$gene %in% CellChatDB$interaction[["receptor"]], TRUE, FALSE)

            # Cell type annotation 
            dea_result <- dplyr::left_join(dea_result, map, by="ident")
            
        }

        i <- as.character(i)
        dea[[i]] <- dea_result
        
    }
    
    # Output
    if(!is.null(file)) {
        
        write.xlsx(dea, paste0(file, ".xlsx"), colNames=TRUE)
        saveRDS(dea, paste0(file, ".rds"))
        
    }

      
    return(dea)

}

##############
### hm_dea ###
##############
hm_dea <- function(dea, so, padj_thr=0.05, log2fc_thr=0.25, top=10, pct_1=NULL, pct_2=NULL, title=NULL, column_name, column_order, row_name=column_name, row_order=column_order, filter_anno=NULL, width=0.25, height=3.0, column_title_rot=0, conserved=FALSE, re_map=FALSE) {
    
    # Filter by annotation 
    if(!is.null(filter_anno)) {dea <- lapply(dea, function(x) {x[x[, filter_anno], ]})}
    
    if(conserved) {
        
        # Compute averages across groups 
        dea <- lapply(dea, function(x) {x$p_val_adj <- x$minimump_p_val; return(x)})
        dea <- lapply(dea, function(x) {x$avg_log2FC <- rowMeans(x[, c("NaCl_avg_log2FC", "CpG_avg_log2FC")]); return(x)})
        dea <- lapply(dea, function(x) {x$pct.1 <- rowMeans(x[, c("NaCl_pct.1", "CpG_pct.1")]); return(x)})
        dea <- lapply(dea, function(x) {x$pct.2 <- rowMeans(x[, c("NaCl_pct.2", "CpG_pct.2")]); return(x)})
        
        # Filter by avg_log2FC sign 
        dea <- lapply(dea, function(x) {x[sign(x$NaCl_avg_log2FC)==sign(x$CpG_avg_log2FC), ]})
        
    }
    
    # Filter by p-value 
    dea <- lapply(dea, function(x) {x[x$p_val_adj <= padj_thr, ]})
    
    # Filter by average log2FC
    dea <- lapply(dea, function(x) {x[x$avg_log2FC>=log2fc_thr, ]})
    
    # Filter by pct.2 
    if(!is.null(pct_2)) {
        
        dea <- lapply(dea, function(x) {x[x$pct.2<=pct_2, ]})
        
    }
    
    # Filter by pct.1 
    if(!is.null(pct_1)) {
        
        dea <- lapply(dea, function(x) {x[x$pct.1>=pct_1, ]})
        
    }
    
    # Order by avg_log2FC
    dea <- lapply(dea, function(x) {x[sign(x$avg_log2FC)==1, ]})
    dea <- lapply(dea, function(x) {x[order(-x$avg_log2FC), ]})

    # Filter by top hits 
    dea <- lapply(dea, function(x) {x[1:top, ]})
    
    # Merge dea results 
    dea <- do.call("rbind", dea)
    dea <- na.omit(dea)
    
    # Re-map cell type names 
    if(re_map) {
        
        dea$cell_type <- NULL
        dea <- dplyr::left_join(dea, re_map_df, by="ident")
        
    }
    
    # Gent normalized count matrix and subset by genes from dea
    mat <- GetAssayData(so, assay="RNA", slot="data", features=dea$gene) %>% as.matrix()
    mat <- mat[rownames(mat) %in% dea$gene, ]
    
    # Scale 
    mat <- t(scale(t(mat)))

    # Order rows by genes
    mat <- mat[dea$gene, ]
    
    # Select column and order from Seurat meta data and order mat
    so@meta.data[[column_name]] <- factor(so@meta.data[[column_name]], levels=column_order)
    column_order_df <- so@meta.data
    column_order_df <- column_order_df[order(column_order_df[[column_name]]), ]
    column_order_df <- column_order_df[column_name]
    
    mat <- mat[, rownames(column_order_df)]
    
    # Breaks 
    breaks <- seq(-quantile(abs(mat), 0.99), quantile(abs(mat), 0.99), length.out=11)
    breaks <- c(breaks[1:4], -0.1, 0, 0.1, breaks[8:11])
    
    color_hm <- c(rev(brewer.pal(11,"RdBu"))[1:5], "#ffffff", rev(brewer.pal(11,"RdBu"))[7:11])
    
    hm <- ComplexHeatmap::Heatmap(
    
        matrix=mat,
        name="z-score", 
        col=colorRamp2(breaks, color_hm), 

        width=ncol(mat)*unit(width, "mm"), 
        height=nrow(mat)*unit(height, "mm"), 

        row_title_gp=gpar(fontsize=12, fontface="plain"),
        column_title_gp=gpar(fontsize=12, fontface="plain"), 

        row_names_gp=grid::gpar(fontsize=8, fontface="plain"), 
        column_names_gp =grid::gpar(fontsize=8, fontface="plain"), 
        
        cluster_rows=TRUE, 
        cluster_columns=TRUE,

        show_row_names=TRUE,
        show_column_names=FALSE, 

        row_gap=unit(1.5, "mm"), 
        column_gap=unit(1.5, "mm"), 

        show_row_dend=FALSE,
        show_column_dend=FALSE,

        row_dend_width=unit(1.5, "mm"),
        column_dend_height=unit(1.5, "mm"),

        row_split=factor(dea[[row_name]], levels=row_order),
        column_split=column_order_df[[column_name]],

        cluster_row_slices=FALSE, 
        cluster_column_slices=FALSE,

        border=TRUE, 

        heatmap_legend_param=list(title_gp=gpar(fontsize=10, fontface="plain"), labels_gp=gpar(fontsize=8), legend_height=unit(10, "mm"), grid_width=unit(2, "mm"))

    ) %>% as.ggplot()
    
    return(hm)

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

########################
### dp_dea_treatment ###
########################
dp_dea_treatment <- function(dea, so, type="receptor", group_by="cell_type_fine", group_by_order=c("cMo (1)", "cMo (2)", "RPMP", "RPM"), padj_thr=0.05, log2fc_thr=0.25, range_min=0, range_max=5) {
    
    require(ggnewscale)
    
    dea <- dea[dea[[type]] & dea$p_val_adj<=padj_thr & abs(dea$avg_log2FC)>=log2fc_thr, ]

    # Extract counts 
    mat <- GetAssayData(so, assay="RNA", slot="data")[dea$gene, ] %>% as.data.frame() %>% add_rownames(var="gene")
    mat <- reshape2::melt(mat, id.vars="gene", value.name="expression", variable.name="cell_id")

    # Combine counts with grouping var and compute variables for plotting 
    mat <- suppressMessages(
        
            dplyr::left_join(mat, so@meta.data[c("cell_id", group_by, "treatment")], by="cell_id") %>%    
            dplyr::group_by_at(c(group_by, "treatment", "gene")) %>% 
            dplyr::summarise(mean_expression=mean(expression, na.rm=TRUE), ratio=sum(expression>0)/n(), cell_count=n()) %>%
            dplyr::ungroup() %>% 
            dplyr::mutate(mean_expression=ifelse(is.nan(mean_expression), 0, mean_expression), ratio=ifelse(ratio==0, NA, ratio))
    
    )

    # Gene group from DEA 
    gene_group <- dplyr::select(dea, gene, avg_log2FC) %>% 
        dplyr::mutate(group=ifelse(sign(avg_log2FC)==1, "CpG", "NaCl")) %>% 
        dplyr::select(gene, group)

    mat <- dplyr::left_join(mat, gene_group, by="gene")

    # Factor levels  
    mat$treatment <- factor(mat$treatment, levels=c("NaCl", "CpG"))
    mat$group <- factor(mat$group, levels=c("NaCl", "CpG"))
    mat[group_by] <- factor(mat[[group_by]], levels=group_by_order)

    # Fill subset data
    mat_nacl <- subset(mat, treatment=="NaCl")
    mat_cpg_subset <- mat_nacl
    mat_cpg_subset$mean_expression <- NA
    mat_cpg_subset$ratio <- NA
    mat_cpg_subset$treatment <- "CpG"
    mat_nacl <- rbind(mat_nacl, mat_cpg_subset)

    mat_cpg <- subset(mat, treatment=="CpG")
    mat_nacl_subset <- mat_cpg
    mat_nacl_subset$mean_expression <- NA
    mat_nacl_subset$ratio <- NA
    mat_nacl_subset$treatment <- "NaCl"
    mat_cpg <- rbind(mat_cpg, mat_nacl_subset)

    # Plot 
    dp <- ggplot(mat) + 

        geom_point(data=mat_nacl, aes_string(x="gene", y=group_by, size="ratio", fill="mean_expression", group="treatment"), position=position_dodge(width=1), alpha=1, shape=21, stroke=0.2) + 
        scale_fill_gradient2("Treatment", limits=c(0, max(abs(mat$mean_expression))), low=color$treatment["NaCl"], mid="white", high=color$treatment["NaCl"]) +

        new_scale("fill") +

        geom_point(data=mat_cpg, aes_string(x="gene", y=group_by, size="ratio", fill="mean_expression", group="treatment"), position=position_dodge(width=1), alpha=1, shape=21, stroke=0.2) + 
        scale_fill_gradient2("CpG", limits=c(0, max(abs(mat$mean_expression))), low=color$treatment["CpG"], mid="white", high=color$treatment["CpG"]) +
        
        scale_size_continuous(name="Fraction", limits=c(0, 1), breaks=c(0, 0.25, 0.5, 0.75, 1.0), range=c(range_min, range_max)) + 
        facet_grid(~group, scale="free", space="free") + 

        xlab("Gene") + ylab("") + ggtitle(type)
    
    return(dp)

}
##############
### vp_dea ###
##############
vp_dea <- function(dea, log2_thold=1, adjpvalue_thold=0.05, top_label=10, title=NULL, conserved=FALSE) {

    if(conserved) {
        
        dea <- dea %>% 
            dplyr::filter(sign(NaCl_avg_log2FC)==sign(CpG_avg_log2FC)) %>% 
            rowwise() %>% 
            dplyr::mutate(avg_log2FC=mean(NaCl_avg_log2FC, CpG_avg_log2FC)) %>% 
            dplyr::rename(p_val_adj=minimump_p_val) %>% 
            as.data.frame()
    }
    
    # Set rownames to genes
    if("gene" %in% colnames(dea)) {rownames(dea) <- dea$gene}
    
    # Check Inf log2FC 
    if(any(is.infinite(dea$avg_log2FC))) {
        
        print(dea[is.infinite(dea$avg_log2FC), ]$gene)
        dea <- dea[!is.infinite(dea$avg_log2FC), ]
        
    }
    
    # Annotate entries significance by log2_thold and adjpvalue_thold
    dea$p_val_adj <- ifelse(dea$p_val_adj == 0, .Machine$double.xmin, dea$p_val_adj)
    dea$sig <- ifelse(abs(dea$avg_log2FC) >= log2_thold & -log10(dea$p_val_adj) >= -log10(adjpvalue_thold), "s", "ns")
    
    # Set color based on significance and direction of dea e.g. positive and negative 
    dea$color <- ifelse(dea$sig == "s" & dea$avg_log2FC > 0, "s_pos", "ns")
    dea$color <- ifelse(dea$sig == "s" & dea$avg_log2FC < 0, "s_neg", dea$color)
    
    color <- c(RColorBrewer::brewer.pal(8, "Set1")[1], "gray", "black", RColorBrewer::brewer.pal(8, "Set1")[2])
    names(color) <- c("s_pos", "ns", "black", "s_neg")
    
    # Create labels based log2FC and p_val_adj
    dea_pos <- dea[dea$avg_log2FC > 0 & dea$sig == "s", ]
    dea_neg <- dea[dea$avg_log2FC < 0 & dea$sig == "s", ]

    pos_labels_log2FC <- dea_pos[rev(order(dea_pos$avg_log2FC)), ][1:top_label, ] %>% rownames()
    neg_labels_log2FC <- dea_neg[order(dea_neg$avg_log2FC), ][1:top_label, ] %>% rownames()
    
    pos_labels_p_val_adj <- dea_pos[order(dea_pos$p_val_adj), ][1:top_label, ] %>% rownames()
    neg_labels_p_val_adj <- dea_neg[order(dea_neg$p_val_adj), ][1:top_label, ] %>% rownames()
    
    pos_labels <- c(pos_labels_log2FC, pos_labels_p_val_adj)
    neg_labels <- c(neg_labels_log2FC, neg_labels_p_val_adj)
    
    # Set labels 
    dea$label <- ifelse(rownames(dea) %in% c(pos_labels, neg_labels), rownames(dea), NA)
    
    # Plot
    library(ggrepel, quietly=TRUE)    
    volcano_plot <- ggplot(dea, aes(x=avg_log2FC, y=-log10(p_val_adj), colour=dea$color, label=label)) + 
        geom_point(shape=16, size=1.5) + 
        geom_vline(aes(xintercept=log2_thold), linetype="dotted", colour="black") +
        geom_vline(aes(xintercept=-log2_thold), linetype="dotted", colour="black") +
        geom_hline(aes(yintercept=-log10(adjpvalue_thold)), linetype="dotted", colour="black") +
        geom_text_repel(segment.color="black", force=20, force_pull=1, max.overlaps=getOption("ggrepel.max.overlaps", default=100), size=3, alpha=1, guide="none", segment.size=0.1, color='black') + 
        xlim(-max(abs(dea$avg_log2FC)), max(abs(dea$avg_log2FC))) + 
        ylim(0, max(-log10(dea$p_val_adj))+max(-log10(dea$p_val_adj))*0.02) + 
        ggtitle(title) + xlab("average log2FC") + ylab("-log10(adj. p-value)") + 
        scale_colour_manual(values=color) + 

        theme(aspect.ratio=1, legend.position="none")
    
    return(volcano_plot)
    
}

####################
### hm_treatment ###
####################
hm_treatment <- function(so, dea, column_name="treatment", cell_type_fine=NULL, top=NULL, log2fc_thr=0.25, padj_thr=0.05, pct_thr=NULL, cluster_columns=TRUE, gene_filter=NULL, width=0.25, height=3.0) {
    
    # Filter genes 
    dea <- dea[!dea$gene %in% gene_filter, ]
    dea <- dea[dea$gene %in% rownames(so), ]
    
    # Filter by padj and log2fc
    dea_nacl <- dea[dea$p_val_adj <= padj_thr & dea$avg_log2FC <= -log2fc_thr, , drop=FALSE]
    dea_cpg <- dea[dea$p_val_adj <= padj_thr & dea$avg_log2FC >= log2fc_thr, , drop=FALSE]
    
    # Filter by pct 
    if(!is.null(pct_thr)) {
        
        dea_nacl <- dea_nacl[dea_nacl$pct.1 <= pct_thr, ]
        dea_cpg <- dea_cpg[dea_cpg$pct.2 <= pct_thr, ]
        
    }

    # Order by padj and log2fc
    dea_nacl <- dea_nacl[with(dea_nacl, order(p_val_adj, avg_log2FC)), ]
    dea_cpg <- dea_cpg[with(dea_cpg, order(p_val_adj, -avg_log2FC)), ]
    
    # Filter top hits
    if(!is.null(top)) {
        
        dea_nacl <- dea_nacl[1:top, ]
        dea_cpg <- dea_cpg[1:top, ]
        
    } 
    
    # Get genes 
    genes_nacl <- dea_nacl$gene %>% na.omit()
    genes_cpg <- dea_cpg$gene %>% na.omit()
    genes <- c(genes_nacl, genes_cpg)    
    
    # Gent normalized count matrix and subset by genes from dea
    mat <- GetAssayData(so, assay="RNA", slot="data", features=genes) %>% as.matrix()
    mat <- mat[rownames(mat) %in% genes, , drop=FALSE]
    
    # Scale 
    mat <- t(scale(t(mat)))
    
    # Order rows by genes
    if(any(genes %in% rownames(mat))==FALSE) {return(NULL)}
    mat <- mat[genes, , drop=FALSE]
    mat <- na.omit(mat)
    
    # If mat is empty return NULL
    if(dim(mat)[1]==0) {return(NULL)}
    
    # Filter genes if all zero after scaling 
    genes_nacl <- genes_nacl[genes_nacl %in% rownames(mat)]
    genes_cpg <- genes_cpg[genes_cpg %in% rownames(mat)]

    # Select column and order from Seurat meta data and order mat
    top_annotation_df <- so@meta.data[, c("treatment", "cell_type_fine")]
    top_annotation_df$treatment <- factor(top_annotation_df$treatment, levels=c("NaCl", "CpG"))
    top_annotation_df$group <- factor(paste0(top_annotation_df$treatment, top_annotation_df$cell_type_fine), levels=c(paste0("NaCl", cell_type_fine), paste0("CpG", cell_type_fine)))
    top_annotation_df <- top_annotation_df[order(top_annotation_df$group), ]
    
    # Top annotaiton based on top_annotation_df
    top_annotation=HeatmapAnnotation(

        treatment=top_annotation_df$treatment, 
        cell_type_fine=top_annotation_df$cell_type_fine, 
        simple_anno_size=unit(4, "mm"), 
        col=list(treatment=unlist(color[["treatment"]]), cell_type_fine=unlist(color[["cell_type_fine"]])), 
        annotation_legend_param=list(title_gp=gpar(fontsize=10, fontface="plain"), labels_gp=gpar(fontsize=8)), 
        show_annotation_name=FALSE

    )
    
    # Order mat by ordered top annotation 
    mat <- mat[, rownames(top_annotation_df), drop=FALSE]
    
    # Cluster by groups 
    cluster_group_df <- split(top_annotation_df, top_annotation_df$group)
    cluster_group_df <- cluster_group_df[unlist(lapply(cluster_group_df, function(x) nrow(x)>0))]
                                                       
    cluster_group_df <- lapply(cluster_group_df, function(x) {rownames(x[column_order(Heatmap(mat[, colnames(mat) %in% rownames(x), drop=FALSE], cluster_columns=TRUE, cluster_rows=FALSE)), , drop=FALSE])})
    mat <- mat[, unlist(cluster_group_df), drop=FALSE]

    # Row split 
    row_split <- data.frame(treatment=factor(c(rep("NaCl", length(genes_nacl)), rep("CpG", length(genes_cpg))), levels=c("NaCl", "CpG")))

    # Color  
    breaks <- seq(-quantile(abs(mat), 0.99), quantile(abs(mat), 0.99), length.out=11)
    breaks <- c(breaks[1:4], -0.1, 0, 0.1, breaks[8:11])
    color_hm <- c(rev(brewer.pal(11,"RdBu"))[1:5], "#ffffff", rev(brewer.pal(11,"RdBu"))[7:11])
    
    hm <- ComplexHeatmap::Heatmap(

        matrix=mat,
        name="z-score", 
        col=colorRamp2(breaks, color_hm), 

        width=ncol(mat)*unit(width, "mm"), 
        height=nrow(mat)*unit(height, "mm"), 

        row_title_gp=gpar(fontsize=12, fontface="plain"),
        column_title_gp=gpar(fontsize=12, fontface="plain"), 

        row_names_gp=grid::gpar(fontsize=8, fontface="plain"), 
        column_names_gp =grid::gpar(fontsize=8, fontface="plain"), 

        cluster_rows=TRUE, 
        cluster_columns=FALSE,

        show_row_names=TRUE,
        show_column_names=FALSE, 

        row_gap=unit(1.5, "mm"), 
        column_gap=unit(1.5, "mm"), 

        show_row_dend=FALSE,
        show_column_dend=FALSE,

        row_dend_width=unit(1.5, "mm"),
        column_dend_height=unit(1.5, "mm"),

        row_split=row_split, 
        column_split=top_annotation_df[["treatment"]],

        top_annotation=top_annotation, 

        cluster_row_slices=FALSE, 
        cluster_column_slices=FALSE,

        border=TRUE, 
        
        use_raster=TRUE, 

        heatmap_legend_param=list(title_gp=gpar(fontsize=10, fontface="plain"), labels_gp=gpar(fontsize=8), legend_height=unit(10, "mm"), grid_width=unit(4, "mm"))



    ) %>% as.ggplot()

    return(hm)
    
}

