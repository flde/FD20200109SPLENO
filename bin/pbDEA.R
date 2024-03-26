####################
### group_select ###
####################
group_select <- function(so, grouping_var, random_effect_var, random_effect_var_cell_min, random_effect_var_min) {
    
    # Get cells grouping information 
    cells_grouping <- so@meta.data[, c(grouping_var, random_effect_var)]
    
    # Count cells per group 
    cells_grouping_n <- cells_grouping %>% 
    
        # Count and filter random_effect_var by random_effect_var_cell_min
        dplyr::group_by_all() %>% 
        dplyr::summarise(random_effect_var_cell_n=n()) %>% 
        dplyr::mutate(random_effect_var_cell_check=ifelse(random_effect_var_cell_n >= random_effect_var_cell_min, TRUE, FALSE)) %>% 
        dplyr::filter(random_effect_var_cell_check) %>% 
    
        # Count and filter grouping_var by random_effect_var_min
        dplyr::group_by(across(grouping_var)) %>% 
        dplyr::mutate(random_effect_var_n=sum(random_effect_var_cell_check)) %>% 
        dplyr::mutate(random_effect_var_check=ifelse(random_effect_var_n >= random_effect_var_min, TRUE, FALSE)) %>% 
        dplyr::filter(random_effect_var_check) %>% 

        # Ungroup 
        dplyr::ungroup() %>% as.data.frame()
    
    # Set cell id filter
    cells_logical <- paste(so@meta.data[[random_effect_var]], so@meta.data[[grouping_var]]) %in% paste(cells_grouping_n[[random_effect_var]], cells_grouping_n[[grouping_var]])
    
    # Set group filter 
    group_logical <- length(unique(cells_grouping_n[[grouping_var]])) > 1
    
    # Filter by cell id and return Seurat object if any cells are found
    if(any(cells_logical) & group_logical) {
        
        so <- subset(so, cells=colnames(so)[cells_logical])
        
    } else {
        
        so <- NULL
    
    }
    
    return(so)
    
}

######################
### feature_select ###
######################
feature_select <- function(so, cnt_min, cell_min, random_effect_var, feature_group_min) {
    
    # Select Seurat components 
    meta <- so@meta.data
    cnt <- GetAssayData(so, assay="RNA", slot="counts")
    
    # Prepare data for split 
    cnt <- t(as.matrix(cnt))
    cnt <- as.data.frame(cnt)
    
    # Split expression matrix by grouping variable into list
    cnt <- split(cnt, f=paste0(so[[grouping_var, drop=TRUE]], "|", so[[random_effect_var, drop=TRUE]]))
    
    # Check within group expression 
    cnt <- lapply(cnt, function(x) {colSums(x>=cnt_min)>=cell_min})
    cnt <- do.call(rbind, cnt)
    
    # Check across group expression 
    cnt_check <- colSums(cnt) >= feature_group_min
    
    # Get genes 
    genes <- colnames(cnt)[cnt_check]

    # Create Seurat object with gene subset 
    cnt <- GetAssayData(so, assay="RNA", slot="counts")
    cnt <- cnt[genes, , drop=FALSE]

    so <- CreateSeuratObject(counts=cnt, meta=meta)
    
    message(paste(unique(so$cell_type_fine), "final number of genes for testing", nrow(so)))

    return(so)
    
}

#######################
### voomlmfit_group ###
#######################
voomlmfit_group <- function(so, grouping_var, random_effect_var, pseudoreplicate_var, contrast_vec, sample_weights) {
    
    # Get Counts 
    cnt <- GetAssayData(so, assay="RNA", slot="counts")
        
    # Prepare count data for split 
    cnt <- t(as.matrix(cnt))
    cnt <- as.data.frame(cnt)
        
    # Make pseudobulks by suming single cell counts
    cnt <- split(cnt, f=paste0(so[[grouping_var, drop=TRUE]], "|", so[[random_effect_var, drop=TRUE]], "|", so[[pseudoreplicate_var, drop=TRUE]]))
    cnt <- lapply(names(cnt), function(i) {x <- data.frame(counts=colSums(cnt[[i]])); colnames(x) <- i; return(x)})
    cnt <- do.call(cbind, cnt)
        
    # Get grouping variables from cnt matrix
    grouping_var_vec <- sapply(strsplit(colnames(cnt), "\\|"), `[[`, 1) # Used for design matrix
    random_effect_var_vec <- sapply(strsplit(colnames(cnt), "\\|"), `[[`, 2) # Used as blocking variable
    
    # Design matrix 
    grouping_var_vec <- as.character(grouping_var_vec)
    
    design <- model.matrix(~0+grouping_var_vec)
    colnames(design) <- gsub("grouping_var_vec", "", colnames(design))
    
    suppressMessages(
        
        fit <- edgeR::voomLmFit(

            counts=cnt, 
            design=design, 
            block=random_effect_var_vec, 
            sample.weights=sample_weights, 
            var.design=design

        )
        
    )
    
    # Set groups for contrast
    contrast_select <- c(contrast_vec %in% grouping_var_vec)
    
    # Set results list
    results <- list()
    
    if(!contrast_select[1]) {
        
        message("Basline samples are missing")
    
    } else if(sum(contrast_select[2:length(contrast_select)])==0) {
        
        message("No groups present to compare to baseline")
        
    } else {
        
        for(i in which(contrast_select)[which(contrast_select)!=1]) {

            # Contrast fit
            contrasts_vec <- paste0(contrast_vec[i], "-", contrast_vec[1])
            contrasts <- limma::makeContrasts(contrasts=contrasts_vec, levels=colnames(design))
            contrasts_fit <- limma::contrasts.fit(fit, contrasts=contrasts)

            # eBayes fit 
            efit <- limma::eBayes(contrasts_fit)

            # Get result table
            result <- limma::topTable(efit, sort.by="P", n=Inf, p.value=1, lfc=0, coef=1)

            # Convert result to Seurat format 
            colnames(result)[1] <- "avg_log2FC"
            colnames(result)[4] <- "p_value"
            colnames(result)[5] <- "p_val_adj"
            
            results[[contrast_vec[i]]] <- result

        }
        
    }
        
    if(length(results)==1) {results <- results[[1]]}

    return(result)
    
}

########################
### voomlmfit_marker ###
########################
voomlmfit_marker <- function(so, grouping_var, pseudobatch_var, sample_weights=FALSE, cnt_min=3, cell_min=3) {
    
    # Get count data 
    cnt <- GetAssayData(so, assay="RNA", slot="counts")
    
    # Get grouping var vector 
    grouping_var_vec_set_i <- so[[grouping_var, drop=TRUE]] %>% unique()
    grouping_var_vec_i <- make.names(grouping_var_vec_set_i)
    so[[grouping_var]] <- make.names(so[[grouping_var, drop=TRUE]])
    so[[pseudobatch_var]] <- make.names(so[[pseudobatch_var, drop=TRUE]])
    
    # Set result list 
    results <- list()
    
    # Iterate over grouping_var_vec_i 
    for(grouping_var_vec_j in grouping_var_vec_i) {
        
        # Filter for genes in grouping_var_vec_j
        cnt_i <- cnt[, so[[grouping_var, drop=TRUE]]==grouping_var_vec_j]
        cnt_i <- cnt_i[rowSums(cnt_i>=cnt_min)>=cell_min, ]
        genes_i <- rownames(cnt_i)

        # Subset counts 
        cnt_i <- cnt[genes_i, ]

        # Prepare count data for split 
        cnt_i <- t(as.matrix(cnt_i))
        cnt_i <- as.data.frame(cnt_i)

        # Make pseudobulks by suming single cell counts
        cnt_i <- split(cnt_i, f=paste0(so[[grouping_var, drop=TRUE]], "|", so[[pseudobatch_var, drop=TRUE]]))
        cnt_i <- lapply(names(cnt_i), function(i) {x <- data.frame(counts=colSums(cnt_i[[i]])); colnames(x) <- i; return(x)})
        cnt_i <- do.call(cbind, cnt_i)

        # Get grouping variables from cnt matrix
        grouping_var_vec <- sapply(strsplit(colnames(cnt_i), "\\|"), `[[`, 1) # Used for design matrix
        pseudobatch_var_vec <- sapply(strsplit(colnames(cnt_i), "\\|"), `[[`, 2) # Used for design matrix
        
        # Design matrix 
        grouping_var_vec <- as.factor(grouping_var_vec)
        grouping_var_vec <- relevel(grouping_var_vec, ref=grouping_var_vec_j) # relevel so that the first design entry is the reference level 
        pseudobatch_var_vec <- as.character(pseudobatch_var_vec)

        design <- model.matrix(~0+grouping_var_vec)
        colnames(design) <- gsub("grouping_var_vec", "", colnames(design))
        
        # voomLmFit 
        fit <- suppressMessages(
            
            edgeR::voomLmFit(

                counts=cnt_i, 
                design=design, 
                sample.weights=FALSE, 
                var.design=design
            
            )
        
        )
        
        # Contrast fit
        contrasts_vec <- paste0(grouping_var_vec_j, "-(", paste0(grouping_var_vec_i[!grouping_var_vec_i %in% grouping_var_vec_j], collapse="+"), ")/", length(grouping_var_vec_i[!grouping_var_vec_i %in% grouping_var_vec_j]))
        contrasts <- limma::makeContrasts(contrasts=contrasts_vec, levels=grouping_var_vec_i)
        contrasts_fit <- limma::contrasts.fit(fit, contrasts=contrasts)

        # eBayes fit 
        efit <- limma::eBayes(contrasts_fit)

        # Get result table
        result <- limma::topTable(efit, sort.by="P", n=Inf, p.value=1, lfc=0, coef=1)

        # Convert result to Seurat format 
        colnames(result)[1] <- "avg_log2FC"
        colnames(result)[4] <- "p_value"
        colnames(result)[5] <- "p_val_adj"
        
        # Add percentage expression 
        cnt_bool <- cnt[rownames(result), ] > 0
        result$pct_1<-rowSums(cnt_bool[, so[[grouping_var, drop=TRUE]]==grouping_var_vec_j])/sum(so[[grouping_var, drop=TRUE]]==grouping_var_vec_j)
        result$pct_2<-rowSums(cnt_bool[, so[[grouping_var, drop=TRUE]]!=grouping_var_vec_j])/sum(so[[grouping_var, drop=TRUE]]!=grouping_var_vec_j)
        
        results[[grouping_var_vec_j]] <- result
        
    }
    
    # Set names 
    names(results) <- grouping_var_vec_set_i
    
    # Return results
    return(results)
    
}

#################
### hm_marker ###
#################
hm_marker <- function(results, so, top=25, grouping_var, width=0.05, height=1, cluster_rows=FALSE, cluster_columns=FALSE) {
    
    # Get grouping_var order from results
    grouping_var_order <- names(results)
    
    # Get top genes
    genes <- lapply(results, function(x) dplyr::filter(x, sign(avg_log2FC)==1, p_val_adj<=0.01) %>% dplyr::mutate(pct_diff=pct_1-pct_2) %>% dplyr::arrange(-pct_diff, avg_log2FC, p_val_adj) %>% dplyr::slice(1:top) %>% rownames)

    # Get normalized data matrix
    mat <- GetAssayData(so, assay="RNA", slot="counts")
    
    # Filter matrix by genes and combine 
    mat <- lapply(genes, function(x) {mat[x, ]})                
    mat <- do.call(rbind, mat)
    mat <- as.matrix(mat)
    
    # Scale and order by column split
    mat <- t(scale(t(mat), center=TRUE, scale=TRUE))
                    
    # Set order 
    so[[grouping_var]] <- factor(so[[grouping_var, drop=TRUE]], levels=grouping_var_order)
                    
    # Order mat columns 
    mat <- mat[, order(so[[grouping_var]])]

    # Splits 
    column_split <- so[[grouping_var, drop=TRUE]][order(so[[grouping_var]])]
    row_split <- factor(rep(names(results), each=top), levels=names(results))
    
    # Color function 
    color_ramp_mat <- viridis::mako(10)
    breaks_mat <- seq(-1, 1, length.out=10)
    color_function_mat <- circlize::colorRamp2(breaks_mat, color_ramp_mat) 
                    
    hm <- ComplexHeatmap::Heatmap(
    
        matrix=mat, 

        col=color_function_mat, 
        na_col="white", 

        width=ncol(mat)*unit(width, "mm"), 
        height=nrow(mat)*unit(height, "mm"), 

        row_title="Cell type marker marker", 
        row_title_gp=gpar(fontsize=18, fontface="bold"), 

        # column_title=NULL, 
        column_title_rot=90, 
        column_title_gp=gpar(fontsize=18, fontface="bold"), 

        row_names_gp=gpar(fontsize=16, fontface="plain"), 
        column_names_gp=gpar(fontsize=18, fontface="plain"), 

        cluster_rows=cluster_rows, 
        cluster_row_slices=FALSE,
        show_row_dend=FALSE,  
        row_split=row_split, 
        row_gap=unit(2, "mm"),
        show_row_names=FALSE,
        row_names_side="left",

        cluster_columns=cluster_columns,
        cluster_column_slices=FALSE, 
        show_column_dend=FALSE, 
        column_split=column_split,
        column_gap=unit(2, "mm"), 
        show_column_names=FALSE, 

        top_annotation=NULL, 
        right_annotation=NULL, 

        rect_gp=gpar(col=NA, lwd=0, alpha=1), 

        heatmap_legend_param=list(title="z-score", title_gp=gpar(fontsize=18, fontface="plain"), labels_gp=gpar(fontsize=16), legend_height=unit(5*6, "mm"), grid_width=unit(6, "mm")), 

        border=TRUE, 
        border_gp=gpar(col="black", lwd=unit(0.5, "mm")), 

        use_raster=FALSE
    
    )
        
    return(hm)
        
}
##############      
### dea_vp ###
##############
dea_vp <- function(dea, log2fc_thr=1, p_adj_thr=0.05, top_label=10, title=NULL, color_neg=RColorBrewer::brewer.pal(8, "Set1")[1], color_pos=RColorBrewer::brewer.pal(8, "Set1")[2]) {

    # Set rownames to genes
    if("gene" %in% colnames(dea)) {rownames(dea) <- dea$gene}
    
    # Annotate entries significance by log2fc_thr and p_adj_thr
    dea$p_val_adj <- ifelse(dea$p_val_adj == 0, .Machine$double.xmin, dea$p_val_adj)
    dea$sig <- ifelse(abs(dea$avg_log2FC) >= log2fc_thr & -log10(dea$p_val_adj) >= -log10(p_adj_thr), "s", "ns")
    
    # Set color based on significance and direction of dea e.g. positive and negative 
    dea$color <- ifelse(dea$sig == "s" & dea$avg_log2FC > 0, "s_pos", "ns")
    dea$color <- ifelse(dea$sig == "s" & dea$avg_log2FC < 0, "s_neg", dea$color)
    
    color <- c(color_neg, "gray", "black", color_pos)
    names(color) <- c("s_neg", "ns", "black", "s_pos")
    
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
    vp <- ggplot(dea, aes(x=avg_log2FC, y=-log10(p_val_adj), fill=dea$color, label=label), alpha=1) + 
    
        geom_point(size=4, shape=21, color="white") + 
        geom_vline(aes(xintercept=log2fc_thr), linetype="dotted", colour="black") +
        geom_vline(aes(xintercept=-log2fc_thr), linetype="dotted", colour="black") +
        geom_hline(aes(yintercept=-log10(p_adj_thr)), linetype="dotted", colour="black") +
        ggrepel::geom_text_repel(segment.color="black", force=20, force_pull=1, max.overlaps=getOption("ggrepel.max.overlaps", default=100), size=5, alpha=1, guide="none", segment.size=0.1, color="black") + 
        xlim(-max(abs(dea$avg_log2FC)), max(abs(dea$avg_log2FC))) +  
        ylim(0, max(-log10(dea$p_val_adj))+1) + 
        ggtitle(title) + xlab("average log2FC") + ylab("-log10(adj. p-value)") + 
        scale_fill_manual(values=color) + 
    
        guides(
            
            color=guide_legend(order=1, title="Group", size=2, keywidth=0.75, keyheight=0.75), 
            alpha="none"
            
        ) + 
    
        theme(

            legend.position="none", 
            aspect.ratio=1

        )
    
    return(vp)
    
}
                    
####################
### Write result ###
####################
write_results <- function(results, basename, subset_xlsx=FALSE) {
    
    # Set names
    ident <- names(results)
    
    # Annotate results 
    results <- lapply(names(results), function(i) {
        
        # Get results 
        result <- results[[i]]
        
        if(nrow(result)>0) {

            # Add meta data to result DEA
            result$gene <- rownames(result)
            result$ident <- i

            # Annotate transcription factors 
            tf <- read.delim("/research/peer/fdeckert/reference/animaltfdb/Mus_musculus_TF.txt")
            result$transcription_factor <- ifelse(result$gene %in% tf$Symbol, TRUE, FALSE)
            
            

            # Annotate receptor and ligands 
            CellChatDB <- CellChat::CellChatDB.mouse
            result$ligand  <- ifelse(result$gene %in% CellChatDB$interaction[["ligand"]], TRUE, FALSE)
            result$receptor  <- ifelse(result$gene %in% CellChatDB$interaction[["receptor"]], TRUE, FALSE)

        }
        
        return(result)

    }
                 )
    
    
    # Set names
    names(results) <- ident
        
    # Output RDS
    saveRDS(results, paste0(basename, "_dea.rds"))
    
    # Output xlsx
    names(results) <- make.names(names(results))
    write.xlsx(results, paste0(basename, "_dea.xlsx"), colNames=TRUE) 
    
    if(subset_xlsx) {
        
        write.xlsx(lapply(results, function(x) {x[x$transcription_factor, ]}), paste0(basename, "_tf_dea.xlsx"), colNames=TRUE) 
        write.xlsx(lapply(results, function(x) {x[x$ligand, ]}), paste0(basename, "_ligand_dea.xlsx"), colNames=TRUE) 
        write.xlsx(lapply(results, function(x) {x[x$receptor, ]}), paste0(basename, "_receptor_dea.xlsx"), colNames=TRUE)
        
    }
 
    
}