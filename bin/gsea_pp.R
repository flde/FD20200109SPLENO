library_load <- suppressMessages(
    
    suppressWarnings(
        
        list(

            library(fgsea), 
            library(msigdbr), 

            library(dplyr)

        )
    )
)

##########################
### GSEA using msigdbr ###
##########################
gsea <- function(dea_res, category="H", subcategory=NULL, gene_set=NULL) {
    
    # Set mgi symbols
    dea_res$mgi_symbol <- rownames(dea_res)
    
    if(is.null(gene_set)) {

        gene_set <- msigdbr::msigdbr(species="mouse", db_species="MM", category=category, subcategory=subcategory)
        
    }
    
    gene_set <- split(gene_set, x=gene_set$gene_symbol, f=gene_set$gs_name)
    
    # Set gene names 
    dea_res <- na.omit(dea_res)
    
    # Make ranks  
    dea_res$p_val_adj <- ifelse(dea_res$p_val_adj==0, min(dea_res$p_val_adj[dea_res$p_val_adj>0]), dea_res$p_val_adj)
    dea_res$sign_log_adj_p_values <- -log10(dea_res$p_val_adj) * sign(dea_res$avg_log2FC)
    
    ranks <- dea_res$sign_log_adj_p_values
    names(ranks) <- dea_res$mgi_symbol
    ranks <- ranks[order(ranks)]
    ranks <- rev(ranks)
    
    # Retain only pathways that overlap with result lsit
    gene_set_filter <- lapply(gene_set, function(x) {sum(dea_res$mgi_symbol %in% x)>=1})
    gene_set <- gene_set[unlist(gene_set_filter)]

    gsea_res <- fgsea::fgsea(
        
        pathways=gene_set,
        stats=ranks,
        # nperm=100000, 
        minSize=1,
        maxSize=Inf
        
    )
    
    return(gsea_res)
    
}

#########################
### ORA using msigdbr ###
#########################
ora <- function(genes, universe, category="MH", subcategory=NULL, gene_set=NULL, min_size=15, max_size=500) {
    
    if(is.null(gene_set)) {

        gene_set <- msigdbr::msigdbr(species="mouse", db_species="MM", category=category, subcategory=subcategory)
        
    }

    print(gene_set %>% dplyr::group_by(gs_name) %>% dplyr::summarise(n=n()) %>% pull(n) %>% summary())
    
    gene_set <- split(gene_set, x=gene_set$gene_symbol, f=gene_set$gs_name)

    gsea_res <- fgsea::fora(
        
        pathways=gene_set,
        genes=genes,
        universe=universe, 
        minSize=min_size,
        maxSize=max_size
        
    )
    
    return(gsea_res)
    
}

##################################
### GSEA using clusterProfiler ###
##################################
gsea_cp <- function(dea_res, ont="BP", pval_thr=0.05) {

    # Set mgi symbols
    dea_res$mgi_symbol <- rownames(dea_res)

    # Make ranks  
    dea_res$p_val_adj <- ifelse(dea_res$p_val_adj==0, min(dea_res$p_val_adj[dea_res$p_val_adj>0]), dea_res$p_val_adj)
    dea_res$sign_log_adj_p_values <- -log10(dea_res$p_val_adj) * sign(dea_res$avg_log2FC)
    
    ranks <- dea_res$sign_log_adj_p_values
    names(ranks) <- dea_res$mgi_symbol
    ranks <- ranks[order(ranks)]
    ranks <- rev(ranks)
    
    # Convert rank names
    conv <- clusterProfiler::bitr(names(ranks), fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Mm.eg.db::org.Mm.eg.db)
    ranks <- ranks[conv$SYMBOL]
    names(ranks) <- conv$ENTREZID

    gsea_res <- clusterProfiler::gseGO(geneList=ranks, OrgDb=org.Mm.eg.db::org.Mm.eg.db, ont=ont, keyType="ENTREZID", minGSSize=10, maxGSSize=500, pAdjustMethod="BH", pvalueCutoff=pval_thr, verbose=FALSE)

    gsea_res <- as.data.frame(gsea_res)

    return(gsea_res)
    
}

#################################
### ORA using clusterProfiler ###
#################################
ora_cp <- function(gene, universe, ont="BP", pval_thr=0.05, qval_thr=0.2) {

    ora_res <- clusterProfiler::enrichGO(gene=gene, OrgDb=org.Mm.eg.db::org.Mm.eg.db, ont=ont, pvalueCutoff=pval_thr, pAdjustMethod="BH", universe=universe, qvalueCutoff=qval_thr, minGSSize=15, maxGSSize=500, readable=FALSE)
    ora_res <- as.data.frame(ora_res)

    return(ora_res)
    
}

#############################
### GSEA GO term simplify ###
#############################
gsea_cp_reduce <- function(gsea_res, ont="BP", reduce_mtx_thr=0.7) {
    
    sem_data <- GOSemSim::godata('org.Mm.eg.db', ont=ont)
    sim_mat <- rrvgo::calculateSimMatrix(unique(gsea_res$ID), semdata=sem_data, orgdb="org.Mm.eg.db", ont=ont, method="Wang")
    
    scores <- gsea_res %>% group_by(ID) %>% summarise(score=min(p.adjust, na.rm=TRUE)) %>% deframe()
    scores <- -log10(scores)
    
    gs_reduce <- rrvgo::reduceSimMatrix(sim_mat, scores, threshold=reduce_mtx_thr, orgdb="org.Mm.eg.db")

    return(gs_reduce)
        
}
