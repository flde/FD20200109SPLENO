########################
### LR result parser ###
########################
lr_res_parse <- function(lr_res, source_i, target_i, sample_group_i, pval_prob_thr) {

    # Filter by sender and target 
    lr_res <- lr_res %>% dplyr::filter(source %in% source_i, target %in% target_i, sample_group %in% sample_group_i)

    # Filter by LR sig 
    lr_sig <- lr_res %>% dplyr::filter(pval<=pval_prob_thr, sample_group %in% sample_group_i) %>% dplyr::pull(interaction_name)
    lr_res <- lr_res %>% dplyr::filter(interaction_name %in% lr_sig)
    
    # Separate ligand complex and set unique interaction_name 
    lr_res <- lr_res %>% separate_rows(., ligand.symbol, sep = ",\\s*") %>% dplyr::mutate(interaction_name=paste0(ligand.symbol, "_", ligand_complex, "_", receptor_complex)) %>% dplyr::distinct()
    
    # Set col name
    col_name <- data.frame(sample_group_i) %>% tibble::rownames_to_column("col_name") %>% dplyr::rename(sample_group=sample_group_i)
    if(all(grepl(":", col_name[[1]]))) {

        col_name <- col_name %>% separate(col_name, into = c("col_name", "bottom_annotation"), sep = ":")
    
    }
    
    lr_res <- dplyr::left_join(lr_res, col_name, by=join_by(sample_group))

    return(lr_res)

}

#########################
### DEA result parser ###
#########################
dea_res_parse <- function(lr_res_i, source_i, sample_group_i) {

    dea_res <- readRDS(paste0("result/dea/scRNAseq/wilcox/", sample_group_i[1], "_vs_", sample_group_i[2], ".rds"))

    # Select DEA result for source 
    dea_res <- dea_res[[source_i]]

    # Create ligand log2fc mat 
    dea_log2fc_i <- dea_res %>% dplyr::filter(gene %in% lr_res_i$ligand.symbol) %>% dplyr::mutate(col_name=names(sample_group_i[1])) %>% 
        dplyr::select(avg_log2FC, gene, col_name) %>% pivot_wider(., names_from=col_name, values_from=avg_log2FC) %>% tibble::column_to_rownames("gene")

    dea_log2fc_i <- dea_log2fc_i %>% tibble::rownames_to_column("ligand.symbol") %>% dplyr::left_join(lr_res_i %>% dplyr::select(interaction_name, ligand.symbol, ligand_complex) %>% dplyr::distinct(), ., by=join_by(ligand.symbol)) %>% 
        dplyr::mutate(sender=source_i)
    
    # Create ligand adj pval mat 
    dea_adj_pval_i <- dea_res %>% dplyr::filter(gene %in% lr_res_i$ligand.symbol) %>% dplyr::mutate(col_name=names(sample_group_i[1])) %>% 
        dplyr::select(p_val_adj, gene, col_name) %>% pivot_wider(., names_from=col_name, values_from=p_val_adj) %>% tibble::column_to_rownames("gene") 

    dea_adj_pval_i <- dea_adj_pval_i %>% tibble::rownames_to_column("ligand.symbol") %>% dplyr::left_join(lr_res_i %>% dplyr::select(interaction_name, ligand.symbol, ligand_complex) %>% dplyr::distinct(), ., by=join_by(ligand.symbol)) %>% 
        dplyr::mutate(sender=source_i)

    # Fix missing results 
    dea_log2fc_i[is.na(dea_log2fc_i)] <- 0
    dea_adj_pval_i[is.na(dea_adj_pval_i)] <- 1

    # Return list
    l_res_i <- list(dea_log2fc_i, dea_adj_pval_i)

    return(l_res_i)
    
}

########################
### DEA result merge ###
########################
dea_res_merge <- function(dea_res_i) {

    l_log2fc_i <- reduce(lapply(dea_res_i, `[[`, 1), full_join, by=join_by(interaction_name, ligand.symbol, ligand_complex, sender)) %>% dplyr::distinct()
    l_log2fc_i[is.na(l_log2fc_i)] <- 0

    l_adj_pval_i <- reduce(lapply(dea_res_i, `[[`, 2), full_join, by=join_by(interaction_name, ligand.symbol, ligand_complex, sender)) %>% dplyr::distinct()
    l_adj_pval_i[is.na(l_adj_pval_i)] <- 1

    dea_res_i <- list(l_log2fc_i, l_adj_pval_i)

    return(dea_res_i)
    
}