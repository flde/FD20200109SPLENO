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
    lr_res <- lr_res %>% separate_rows(., ligand_symbol, sep = ",\\s*") %>% dplyr::mutate(interaction_name=paste0(ligand_symbol, "_", receptor_symbol)) %>% dplyr::distinct()
    
    # Set col name
    col_name <- data.frame(sample_group_i) %>% tibble::rownames_to_column("col_name") %>% dplyr::rename(sample_group=sample_group_i)
    if(all(grepl(":", col_name[[1]]))) {

        col_name <- col_name %>% separate(col_name, into = c("col_name", "bottom_annotation"), sep = ":")
    
    }
    
    lr_res <- dplyr::left_join(lr_res, col_name, by=join_by(sample_group))

    return(lr_res)

}

########################
### LR result strech ###
########################
lr_res_stretch <- function(lr_res) {

    lr_res <- lr_res %>% tidyr::separate_rows(ligand_symbol, sep="_") %>% tidyr::separate_rows(receptor_symbol, sep="_")

    return(lr_res)
    
}