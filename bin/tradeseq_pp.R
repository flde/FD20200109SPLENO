library_load <- suppressMessages(
    
    list(
            library(Seurat), 
            library(SeuratWrappers), 
            library(tradeSeq), 
            library(SingleCellExperiment), 

            
            # IHW 
            library(IHW),
            
            # Reactome GSEA/ORA and data base 
            library(ReactomePA),
            library(clusterProfiler), # Symbol to ID
            library(org.Mm.eg.db),

            # Data 
            library(dplyr), 
            library(tidyverse), 
            library(data.table)
        
    )
    
)

##############################
### PT order rows of a mat ###
##############################
pt_order <- function(mat) {

    row_split <- apply(mat, 1, which.max)
    row_max <- apply(mat, 1, max)
    
    row_split_order <- data.frame(row_names=names(row_split), row_split=row_split, row_max=row_max) %>% dplyr::arrange(row_split)
    row_split_order <- split(row_split_order, row_split_order$row_split)
    
    row_split_order <- lapply(seq_along(row_split_order), function(i) {
        
        if(i<length(row_split_order)) {
        
            row_split_order_i <- row_split_order[[i]]
            mat_i <- mat[rownames(row_split_order_i), (i+1), drop=FALSE]
            colnames(mat_i) <- "row_max_i"
        
            row_split_order_i <- cbind(row_split_order_i, mat_i) %>% dplyr::arrange(row_max_i, desc(row_max))
                
        } else if(i>1 & i<length(row_split_order)) {
        
            row_split_order_i <- row_split_order[[i]]
            mat_i <- mat[rownames(row_split_order_i), (i+1), drop=FALSE]
            mat_j <- mat[rownames(row_split_order_j), (i-1), drop=FALSE]
            colnames(mat_j) <- "row_max_j"
        
            row_split_order_i <- cbind(row_split_order_i, mat_i, mat_j) %>% dplyr::arrange(desc(row_max_j), row_max, row_max_i)
                
        } else {
        
            row_split_order_i <- row_split_order[[i]]
            mat_i <- mat[rownames(row_split_order_i), (i-1), drop=FALSE]
            colnames(mat_i) <- "row_max_i"
        
            row_split_order_i <- cbind(row_split_order_i, mat_i) %>% dplyr::arrange(desc(row_max_i), row_max)
                
        }
            
    }
                             )

    row_split_order <- do.call(rbind, row_split_order)
    row_split_order <- setNames(row_split_order$row_split, row_split_order$row_names)

    return(row_split_order)
    
}

######################
### PT split genes ###
######################
pt_split <- function(fitgam, genes, condition_vec, n_points=101) {

    # Get knots positions to test 
    pt_knots <- readRDS("result/lineage/pt_knots.rds")[["pt_knots"]]

    # Get smoothed data 
    mat <- predictSmooth(fitgam, gene=genes, nPoints=n_points, tidy=TRUE)

    # filter by time point bins 
    mat <- mat %>% dplyr::filter(time %in% pt_knots)
    
    if(length(condition_vec)==1) {
        
        mat <- mat %>% dplyr::filter(condition %in% condition_vec)
        mat <- mat %>% dplyr::mutate(time=paste0(time, ":", condition)) %>% dplyr::select(-lineage, -condition) %>% pivot_wider(., names_from=time, values_from=yhat) %>% tibble::column_to_rownames("gene")

        mat <- t(apply(as.matrix(mat), 1, scales::rescale))
        
    } else {

        mat <- lapply(condition_vec, function(i) {

            mat <- mat %>% dplyr::filter(condition %in% i)
            mat <- mat %>% dplyr::mutate(time=paste0(time, ":", condition)) %>% dplyr::select(-lineage, -condition) %>% pivot_wider(., names_from=time, values_from=yhat) %>% tibble::column_to_rownames("gene")
            
            mat <- t(apply(as.matrix(mat), 1, scales::rescale))
    
            
        }
                     )

        mat <- Reduce("+", mat) / length(mat)
        
    }

    # Set split categories
    mat <- list(
        
        MEP=rowMeans(mat[, c(1, 2)], na.rm=TRUE), 
        ProEB=mat[, c(3)], 
        EB=rowMeans(mat[, c(4, 5)], na.rm=TRUE)
        
    )
    
    mat <- do.call(cbind, mat)

    # Order genes by splits
    row_split_order <- pt_order(mat)
    
    return(row_split_order)
    
}

######################
### PT split genes ###
######################
pt_diff_split <- function(fitgam, genes, condition_qry, condition_ref, n_points=101) {

    # Get knots positions to test 
    pt_knots <- readRDS("result/lineage/pt_knots.rds")[["pt_knots"]]

    # Get smoothed data 
    mat <- predictSmooth(fitgam, gene=genes, nPoints=n_points, tidy=TRUE)

    # Get pt knots 
    mat <- mat %>% dplyr::filter(time %in% pt_knots)

    # Average of query scaled matrix
    mat <- mat %>% dplyr::filter(time %in% pt_knots)
    
    mat_qry <- lapply(condition_qry, function(i) {

        mat <- mat %>% dplyr::filter(condition %in% i)
        mat <- mat %>% dplyr::mutate(time=paste0(time, ":", condition)) %>% dplyr::select(-lineage, -condition) %>% pivot_wider(., names_from=time, values_from=yhat) %>% tibble::column_to_rownames("gene")
    
            
    }
                 )

    mat_qry <- Reduce("+", mat_qry) / length(mat_qry)
    
    # Average of reference scaled matrix
    mat_ref <- lapply(condition_ref, function(i) {

        mat <- mat %>% dplyr::filter(condition %in% i)
        mat <- mat %>% dplyr::mutate(time=paste0(time, ":", condition)) %>% dplyr::select(-lineage, -condition) %>% pivot_wider(., names_from=time, values_from=yhat) %>% tibble::column_to_rownames("gene")
    
            
    }
                 )

    mat_ref <- Reduce("+", mat_ref) / length(mat_ref)

    # log2FC matrix
    mat_fc <- log2((mat_qry+1)/(mat_ref+1))

    # Set neg FC to zero to not influence the order of genes
    # mat_fc[mat_fc<0] <- 0
    
    # Set split categories
    mat_fc <- list(
        
        MEP=rowMeans(mat_fc[, c(1, 2)], na.rm=TRUE), 
        ProEB=mat_fc[, c(3)], 
        EB=rowMeans(mat_fc[, c(4, 5)], na.rm=TRUE)
        
    )

    mat_fc <- do.call(cbind, mat_fc)

    # Scale to focus on peak FC
    mat_fc <- t(apply(as.matrix(mat_fc), 1, scales::rescale))

    # Order genes by splits
    row_split_order <- pt_order(mat_fc)
    
    return(row_split_order)
    
}

#################################################
### Prepare TF TG matrix for heatmap plotting ###
#################################################
tf_tg_pp <- function(tf_tg_mat) {

    tf_tg_mat <- tf_tg_mat %>% mutate(value=1) %>% pivot_wider(names_from=TF, values_from=value, values_fill=0) %>% tibble::column_to_rownames("target_genes")


    return(tf_tg_mat)

}
##############################################
### Prune TF TG matrix by split regulation ###
##############################################
tf_tg_prune <- function(tf_tg_mat, split, tf_ratio=0, tf_sum=0, tf_fill=TRUE) {

    tf_prune <- lapply(unique(split), function(i) {

        genes_i <- split[split==i] %>% names()
        tf_tg_mat_i <- tf_tg_mat[rownames(tf_tg_mat) %in% genes_i, ]
        
        tf_tg_mat_i <- tf_tg_mat_i[, (colSums(tf_tg_mat_i)/nrow(tf_tg_mat_i))>=tf_ratio]
        tf_tg_mat_i <- tf_tg_mat_i[, colSums(tf_tg_mat_i)>=tf_sum]

        tf_prune_i <- colnames(tf_tg_mat_i)

        return(tf_prune_i)
        
    }
          )

    tf_prune <- do.call(c, tf_prune) %>% unique()

    tf_tg_mat <- tf_tg_mat[, tf_prune]

    # Fill missing rows to match split target genes 
    if(tf_fill) {

        genes <- names(split)
        
        tf_tg_mat <- rbind(
    
            tf_tg_mat, 
            matrix(0, nrow=length(genes[!genes %in% rownames(tf_tg_mat)]), ncol=ncol(tf_tg_mat), dimnames=list(genes[!genes %in% rownames(tf_tg_mat)], colnames(tf_tg_mat)))
        
        )
    
    }

    return(tf_tg_mat)

}

#######################
### TF column order ###
#######################
tf_tg_order <- function(tf_tg_mat, split) {

    split <- split[names(split) %in% rownames(tf_tg_mat)]
    split_mat <- tf_tg_mat[names(split), ]
    split_mat <- split(split_mat, split)
    split_mat <- lapply(split_mat, colSums) %>% do.call(rbind, .)
    
    tf_tg_order <- pt_order(t(split_mat))
    
    tf_tg_mat <- tf_tg_mat[, names(tf_tg_order)]

    return(tf_tg_mat)
    
}

#############################################################
### Select top hits by smoothed log2FC in sliding windows ###
#############################################################
pt_diff_smooth <- function(fitgam, genes, condition_qry, condition_ref, n_points=21) {

    # Get knots positions to test 
    pt_knots <- seq(0, 1, length.out=n_points)
    pt_knots <- c(0, pt_knots[seq(2, length(pt_knots), by=2)], 1)
    
    # Get genes with positive log2FC
    mat <- predictSmooth(fitgam, genes, nPoints=n_points, tidy=TRUE) %>% 

        # Subset by knots that match wilcox time bin center 
        dplyr::filter(time %in% pt_knots) %>%

        # Filter for condition 
        dplyr::filter(condition %in% c(condition_qry, condition_ref)) %>% 

        dplyr::group_by(gene) %>% dplyr::mutate(yhat=scales::rescale(yhat)) %>% dplyr::ungroup() %>%

        # Set reference condition 
        dplyr::mutate(condition=ifelse(condition==condition_qry, "condition_1", "condition_2")) %>%
    
        # Transform smooth expression to matrix
        pivot_wider(names_from=condition, values_from=yhat) %>% 
        dplyr::group_by(gene) %>% 
    
        # Sliding window 
        dplyr::arrange(gene, time) %>% 
        dplyr::mutate(time_idx=row_number()) %>% 
        dplyr::mutate(diff_smooth=condition_1-condition_2) %>% 
        dplyr::rename(pt_bin=time) %>% dplyr::mutate(pt_bin=factor(pt_bin))

    mat$condition_qry <- condition_qry
    mat$condition_ref <- condition_ref

    return(mat)

}
#########################
### Presto Wilcox test###
#########################
wilcox <- function(so, ident, ident_1=NULL, ident_2=NULL, only_pos=FALSE, avg_log2FC_threshold=0, min_pct=0, test_use="wilcox") {
    
    # Drop empty levels 
    so@meta.data <- droplevels(so@meta.data)
    
    # Check number of cells 
    n_cells_1 <- sum(so@meta.data[[ident]]==ident_1)
    n_cells_2 <- sum(so@meta.data[[ident]]==ident_2)

    check_1 <- n_cells_1 >= 3
    check_2 <- n_cells_2 >= 3

    if(check_1 & check_2) {
        
        so <- suppressMessages(NormalizeData(so))
        so <- SetIdent(so, value=ident)
        res <- RunPresto(so, ident.1=ident_1, ident.2=ident_2, logfc.threshold=avg_log2FC_threshold, min.pct=min_pct, only.pos=only_pos, test.use=test_use)

        # Adjusted p-value with IHW 
        res$mean_exp <- rowMeans(GetAssayData(so, assay="RNA", slot="data")[rownames(so), ])
        # res$p_val_adj <- IHW::adj_pvalues(IHW::ihw(res$p_val ~ res$mean_exp, alpha=0.05))

        # Annotate results 
        res$gene <- rownames(res)

        # N cells per group 
        res$n_cells_1 <- n_cells_1
        res$n_cells_2 <- n_cells_2

        return(res)
        
    } else {
        
        return(NULL)
        
    }

}

##############################################
### Select top hits by Wilcox in time bins ###
##############################################
pt_diff_wilcox <- function(fitgam, genes, condition_qry, condition_ref) {

    # Get dp bins for testing 
    pt_bins <- seq(0, 1, length.out=11)

    # Get pt knots for labeling 
    pt_knots <- seq(0, 1, length.out=21)
    pt_knots <- pt_knots[seq(2, length(pt_knots), by=2)]

    # Get counts
    counts <- assay(fitgam)[genes, , drop=FALSE]
    
    # Get meta data
    meta <- data.frame(condition=fitgam$tradeSeq$conditions, dpt=fitgam$crv$pseudotime)
    rownames(meta) <- colnames(fitgam)
    
    # Create Seurat object
    so <- CreateSeuratObject(counts=counts, meta.data=meta, assay="RNA")
    
    # Get pt bins 
    so$pt_bin <- cut(so$dpt, breaks=pt_bins, labels=pt_knots, include.lowest=TRUE) 
    
    # Split by bins
    so <- SplitObject(so, split.by="pt_bin")

    so <- lapply(so, function(x) {
        
        # Split by condition
        x_1 <- x[, x$condition == condition_qry]
        x_2 <- x[, x$condition == condition_ref]
        
        # Get minimal cell count
        min_n <- min(ncol(x_1), ncol(x_2))
        
        # Subsample each
        set.seed(42)
        x_1_sub <- x_1[, sample(colnames(x_1), min_n)]
        x_2_sub <- x_2[, sample(colnames(x_2), min_n)]
        
        # Merge back
        x <- merge(x_1_sub, x_2_sub)
        x <- JoinLayers(x)

        return(x)
    }
                )

    # Run Wilcox 
    res_dea <- lapply(so, function(so) {
        
        dea_res_i <- wilcox(
            
            so=so, 
            ident="condition", 
            ident_1=condition_qry, 
            ident_2=condition_ref, 
            only_pos=FALSE, 
            avg_log2FC_threshold=0, 
            min_pct=0, 
            test_use="wilcox"
        
        )

        dea_res_i$pt_bin <- so$pt_bin[1]
        dea_res_i$gene <- rownames(dea_res_i)

        dea_res_i <- dea_res_i %>% dplyr::rename(log2FC_wilcox=avg_log2FC, p_val_adj_wilcox=p_val_adj)

        return(dea_res_i)
    
    }
                     )

    res_dea <- do.call(rbind, res_dea)

    return(res_dea)
    
}

#################################################
### Prepare TF TG matrix for heatmap plotting ###
#################################################
enrichment_pp <- function(enrichment_res, pvalue_thr=0.1) {

    if("geneID" %in% colnames(enrichment_res)) {

        enrichment_res <- enrichment_res %>% 
            dplyr::rename("gene"="geneID") %>%
            dplyr::filter(pvalue<=pvalue_thr) %>% 
            dplyr::select(top_level, ID, Description, gene) %>% 
            separate_rows(., gene, sep="/\\s*") %>% 
            dplyr::rename(pathway=Description, gene=gene) %>% 
            dplyr::distinct()
        
    }

    if("core_enrichment" %in% colnames(enrichment_res)) {
        
        # Convert ora results to mat
        enrichment_res <- enrichment_res %>%
            dplyr::rename("gene"="core_enrichment") %>%
            dplyr::filter(pvalue<=pvalue_thr) %>% 
            dplyr::select(top_level, ID, Description, gene) %>% 
            separate_rows(., gene, sep="/\\s*") %>% 
            dplyr::mutate(gene=id_to_symbol(gene)) %>% 
            dplyr::rename(pathway=Description, gene=gene) %>% 
            dplyr::distinct()
        
    }
    
    if(boost) {

        library(reactome.db)
        reactome_id_genes <- as.list(reactome.db::reactomePATHID2EXTID)
    
        reactome_boost <- reactome_id_genes[enrichment_res$ID] %>% stack %>% dplyr::rename(ENTREZID=values, ID=ind) %>% dplyr::mutate(gene=id_to_symbol(ENTREZID)) %>% dplyr::select(ID, gene) %>% dplyr::distinct()
        
        reactome_boost <- reactome_boost[reactome_boost$gene %in% diff_genes, ] # Only take boosting genes that are perturbed 
        reactome_boost <- dplyr::left_join(dplyr::select(enrichment_res, -gene), reactome_boost, by=join_by(ID)) %>% dplyr::distinct()

        enrichment_res <- rbind(enrichment_res, reactome_boost) %>% na.omit() %>% dplyr::distinct() 
        
    }
    
    enrichment_mat <- enrichment_res %>% dplyr::select(top_level, gene) %>% dplyr::distinct() %>% mutate(value=1) %>% pivot_wider(names_from=top_level, values_from=value, values_fill=0) %>% 
        tibble::column_to_rownames("gene")
    
}
    
################
### pt_ranks ###
################
pt_rank <- function(fitgam, association_res, ptpg_res, condition_ref, n_points, log2fc_thr, log2fc_sw, exp_min, exp_min_sw, exp_max) {

    # Get pt diff results 
    pt_diff_res <- pt_diff(fitgam, association_res, ptpg_res$gene, condition_ref=condition_ref, n_points=n_points, log2fc_thr=log2fc_thr, log2fc_sw=log2fc_sw, exp_min=exp_min, exp_min_sw=exp_min_sw, exp_max=exp_max)

    # Add log2fc to ptpg results 
    ptpg_res <- ptpg_res %>% dplyr::mutate(log2fc=ifelse(gene %in% pt_diff_res$pt_diff_1, +1, ifelse(gene %in% pt_diff_res$pt_diff_2, -1, 0))) # CHECK pt_diff_res$pt_diff_1

    # Set undecided pvalues to 0 
    ptpg_res <- ptpg_res[abs(ptpg_res$log2fc)>0, ]

    # Make ranks 
    ptpg_res$pval <- ifelse(ptpg_res$pval==0, min(ptpg_res$pval[ptpg_res$pval>0]), ptpg_res$pval)
    ptpg_res$sign_log_adj_p_values <- -log10(ptpg_res$pval) * sign(ptpg_res$log2fc)
    
    ranks <- ptpg_res$sign_log_adj_p_values
    names(ranks) <- ptpg_res$gene
    ranks <- ranks[order(ranks)]
    ranks <- rev(ranks)

    return(ranks)
    
}

############################
### Prepare Reactome mat ###
############################
pt_reactome_pp <- function(mat, adj_pval_thr, boost=FALSE, diff_genes=NULL) {

    # Convert GSEA results to mat
    if("geneID" %in% colnames(mat)) {
        
        mat <- mat %>% 
            dplyr::rename("gene"="geneID") %>%
            dplyr::filter(pvalue<=adj_pval_thr) %>% # Fix filter
            dplyr::select(top_level, ID, Description, gene) %>% 
            separate_rows(., gene, sep="/\\s*") %>% 
            dplyr::rename(pathway=Description, gene=gene) %>% 
            dplyr::distinct()
            
    }

    # Convert ORA results to mat
    if("core_enrichment" %in% colnames(mat)) {
            
        mat <- mat %>%
            dplyr::rename("gene"="core_enrichment") %>%
            dplyr::filter(pvalue<=0.1) %>% 
            dplyr::select(top_level, ID, Description, gene) %>% 
            separate_rows(., gene, sep="/\\s*") %>% 
            dplyr::mutate(gene=id_to_symbol(gene)) %>% 
            dplyr::rename(pathway=Description, gene=gene) %>% 
            dplyr::distinct()
            
        }

    # Boost pathway genes
    if(boost & !(is.null(diff_genes))) {
    
        library(reactome.db)
        reactome_id_genes <- as.list(reactome.db::reactomePATHID2EXTID)
        
        reactome_boost <- reactome_id_genes[mat$ID] %>% stack %>% dplyr::rename(ENTREZID=values, ID=ind) %>% dplyr::mutate(gene=id_to_symbol(ENTREZID)) %>% dplyr::select(ID, gene) %>% dplyr::distinct()
            
        reactome_boost <- reactome_boost[reactome_boost$gene %in% diff_genes, ] # Only take boosting genes that are perturbed 
        reactome_boost <- dplyr::left_join(dplyr::select(mat, -gene), reactome_boost, by=join_by(ID)) %>% dplyr::distinct()
    
        mat <- rbind(mat, reactome_boost) %>% na.omit() %>% dplyr::distinct() 
            
    }

    return(mat)

}



























