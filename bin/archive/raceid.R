#######################
### RaceID Pipe #######
#######################
raceid_pp <- function(so, suffix, CGenes=NULL, compute=TRUE) {
    
        if(compute) {
        
            # Create SCseq class from Seurat object 
            cnt <- GetAssayData(so, assay="RNA", slot="counts")
            cnt <- cnt[rowSums(cnt)>0, ]
            sc <- SCseq(expdata=cnt)

            # Filter SCseq object  
            sc <- filterdata(sc, mintotal=1)

            expData <- getExpData(sc)
            
            # Cell cycle score for regVar 
            s_score <- so$msS_RNA
            g2m_score <- so$msG2M_RNA
            
            # regVar 
            regVar <- data.frame(s_score=s_score, g2m_score=g2m_score)
            rownames(regVar) <- colnames(expData)
            
            # Batch 
            batch <- so$sample_group 
            names(batch) <- colnames(so)
            
            # pruneKnn
            res <- pruneKnn(expData, no_cores=8, batch=batch, regVar=regVar, alpha=NULL, bmethod="harmony")

            reticulate::use_python("/home/fdeckert/bin/miniconda3/envs/p.3.8.12-FD20200109SPLENO/bin/python3", required=TRUE)
            cl <- RaceID::graphCluster(res, pvalue=0.01, use.leiden=TRUE, leiden.resolution=1)

            sc <- updateSC(sc, res=res, cl=cl)
            sc <- compumap(sc, min_dist=0.3)
            probs <- transitionProbs(res, cl, pvalue=0.01)

            raceid <- list(sc=sc, expData=expData, res=res, probs=probs)

            saveRDS(raceid, paste0("data/object/raceid/raceid", suffix, ".rds"))

            return(raceid)
        
        }

    else {
        
        raceid <- readRDS(paste0("data/object/raceid/raceid", suffix, ".rds"))
        
        return(raceid)
    
    }
    
}

#################################
### RaceID Pipe to Seurat #######
#################################
raceid_to_seurat <- function(so, raceid){
    
    sc <- raceid[["sc"]]
    
    umap_varid <- do.call("cbind", sc@umap@.Data)
    umap_varid <- as.matrix(umap_varid)
    colnames(umap_varid) <- c("UMAP_1", "UMAP_2")
    rownames(umap_varid) <- names(sc@cpart)
    so@reductions[["umap_varid"]] <- CreateDimReducObject(embeddings= umap_varid, key="UMAPvarid_", assay="RNA")
    
    umap_varid <- do.call("cbind", sc@fr@.Data)
    umap_varid <- as.matrix(umap_varid)
    colnames(umap_varid) <- c("FR_1", "FR_2")
    rownames(umap_varid) <- names(sc@cpart)
    so@reductions[["fr_varid"]] <- CreateDimReducObject(embeddings= umap_varid, key="FRvarid_", assay="RNA")
    
    varid <- data.frame(varid_clusters=sc@cluster$kpart)
    rownames(varid) <- names(sc@cpart)
    
    if("varid_clusters" %in% colnames(so@meta.data)) {so$varid_clusters <- NULL}

    so <- AddMetaData(so, varid)
    
    return(so)
    
}