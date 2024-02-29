feature_select <- function(so, cnt_min, cell_min, grouping_var, random_effect_var, feature_group_min) {
    
    # Select Seurat components 
    meta <- so@meta.data
    cnt <- GetAssayData(so, assay="RNA", slot="counts")
    
    # Prepare data for split 
    cnt <- t(as.matrix(cnt))
    cnt <- as.data.frame(cnt)
    
    # Split expression matrix by grouping variable into list
    cnt <- split(cnt, f=paste(so[[grouping_var, drop=TRUE]], so[[random_effect_var, drop=TRUE]]))
    
    # Check within group expression 
    cnt <- lapply(cnt, function(x) {colSums(x>=cnt_min)>cell_min})
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

voom_lm_fit <- function(so, grouping_var, random_effect_var, sample_weights) {
    
    # Get Counts 
    cnt <- GetAssayData(so, assay="RNA", slot="counts")
        
    # Prepare count data for split 
    cnt <- t(as.matrix(cnt))
    cnt <- as.data.frame(cnt)
        
    # Make pseudobulks by suming single cell counts
    cnt <- split(cnt, f=paste0(so[[grouping_var, drop=TRUE]], "|", so[[random_effect_var, drop=TRUE]]))
    cnt <- lapply(names(cnt), function(i) {x <- data.frame(counts=colSums(cnt[[i]])); colnames(x) <- i; return(x)})
    cnt <- do.call(cbind, cnt)
        
    # Get grouping variables from cnt matrix
    grouping_var_vec <- sapply(strsplit(colnames(cnt), "\\|"), `[[`, 1) # Used for design matrix
    random_effect_var_vec <- sapply(strsplit(colnames(cnt), "\\|"), `[[`, 2) # Used as blocking variable
    
    print(grouping_var_vec)
    print(random_effect_var_vec)
 
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
    
    return(fit)
    
}

ebayes <- function(fit) {
    
    # Get grouping var and design from fit
    grouping_var_vec <- colnames(fit$design)

    design <- model.matrix(~0+grouping_var_vec)
    colnames(design) <- gsub("grouping_var_vec", "", colnames(design))
    
    # Define result list
    result <- list()
    
    # Set groups for contrast
    grouping_var_vec_1 <- c("NaCl", "CpG")

    grouping_select <- c(grouping_var_vec_1 %in% grouping_var_vec)
    
    if(!grouping_select[1]) {
        
        message("Basline samples are missing")
    
    } else if(sum(grouping_select[2:length(grouping_select)])==0) {
        
        message("No groups present to compare to baseline")
        
    } else {
        
        for(i in which(grouping_select)[which(grouping_select)!=1]) {

            # Contrast fit
            contrasts_vec <- paste0(grouping_var_vec_1[i], "-", grouping_var_vec_1[1])
            contrasts <- limma::makeContrasts(contrasts=contrasts_vec, levels=colnames(design))
            contrasts_fit <- limma::contrasts.fit(fit, contrasts=contrasts)

            # eBayes fit 
            efit <- limma::eBayes(contrasts_fit)

            # Get result table
            result_df <- limma::topTable(efit, sort.by="P", n=Inf, p.value=1, lfc=0, coef=1)

            # Convert results to Seurat format 
            colnames(result_df)[1] <- "avg_log2FC"
            colnames(result_df)[4] <- "p_value"
            colnames(result_df)[5] <- "p_val_adj"

            result[[contrasts_vec]] <- result_df

        }
        
    }
        
    # 
    if(length(grouping_var_vec_1) < 3) {result <- result[[1]]}

    return(result)
    
}

write_dea <- function(dea, file) {
    
    # Set names
    ident <- names(dea)
    
    # Annotate results 
    dea <- lapply(names(dea), function(i) {
        
        # Get results 
        dea_result <- dea[[i]]
        
        if(nrow(dea_result)>0) {

            # Add meta data to result DEA
            dea_result$gene <- rownames(dea_result)
            dea_result$ident <- i

            # Annotate transcription factors 
            tf <- read.delim("/research/peer/fdeckert/reference/animaltfdb/Mus_musculus_TF.txt")
            dea_result$transcription_factor <- ifelse(dea_result$gene %in% tf$Symbol, TRUE, FALSE)
            
            

            # Annotate receptor and ligands 
            CellChatDB <- CellChat::CellChatDB.mouse
            dea_result$ligand  <- ifelse(dea_result$gene %in% CellChatDB$interaction[["ligand"]], TRUE, FALSE)
            dea_result$receptor  <- ifelse(dea_result$gene %in% CellChatDB$interaction[["receptor"]], TRUE, FALSE)

        }
        
        return(dea_result)

    }
                 )
    
    
    # Set names
    names(dea) <- ident
        
    # Output RDS
    saveRDS(dea, paste0(file, "_dea.rds"))
    
    # Output xlsx
    names(dea) <- gsub("\\.", "", make.names(names(dea)))
    write.xlsx(dea, paste0(file, "_dea.xlsx"), colNames=TRUE) 
    write.xlsx(lapply(dea, function(x) {x[x$transcription_factor, ]}), paste0(file, "_tf_dea.xlsx"), colNames=TRUE) 
    write.xlsx(lapply(dea, function(x) {x[x$ligand, ]}), paste0(file, "_ligand_dea.xlsx"), colNames=TRUE) 
    write.xlsx(lapply(dea, function(x) {x[x$receptor, ]}), paste0(file, "_receptor_dea.xlsx"), colNames=TRUE) 
    
}