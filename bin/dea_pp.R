library_load <- suppressMessages(
    
    suppressWarnings(
        
        list(
        
            library(dplyr)

        )
    )
)


feature_select <- function(so, pct_min=0, cnt_min=3, cell_min=1) {

    # Get count matrix 
    cnt <- GetAssayData(so, assay="RNA", layer="counts")

    # Filter by percentage 
    cnt <- cnt[(100*rowSums(cnt>0)/ncol(cnt))>=pct_min, ]

    # Filer by expression 
    cnt <- cnt[rowSums(cnt>=cnt_min)>=cell_min, ]

    # return genes 
    genes <- rownames(cnt)
        
    return(genes)
    
}

wilcox <- function(so, ident, ident_1=NULL, ident_2=NULL, only_pos=FALSE, avg_log2FC_threshold=0, pct_min=0, cnt_min=0, cell_min=0, test_use="wilcox", verbose=TRUE) {

    # Select genes per group 
    genes_ident_1 <- feature_select(so[, so$group==ident_1], pct_min=pct_min, cnt_min=cnt_min, cell_min=cell_min)
    genes_ident_2 <- feature_select(so[, so$group==ident_2], pct_min=pct_min, cnt_min=cnt_min, cell_min=cell_min)

    genes <- union(genes_ident_1, genes_ident_2)

    genes <- genes[!genes %in% c("Igha", "Igkc", "Tcrg-C1")]

    if(verbose) message(paste0(so$celltype_low[1], " genes selected: ", length(genes)))

    so <- so[genes, ]
    
    # Drop empty levels 
    so@meta.data <- droplevels(so@meta.data)
    
    # Check number of cells 
    n_cells_1 <- sum(so@meta.data[[ident]]==ident_1)
    n_cells_2 <- sum(so@meta.data[[ident]]==ident_2)

    check_1 <- n_cells_1 >= 3
    check_2 <- n_cells_2 >= 3
    
    if(check_1 & check_2) {
        
        so <- SetIdent(so, value=ident)
        so <- NormalizeData(so)
        res <- RunPresto(so, ident.1=ident_1, ident.2=ident_2, logfc.threshold=avg_log2FC_threshold, min.pct=0, only.pos=only_pos, test.use=test_use)

        # Adjusted p-value with IHW 
        res$mean_exp <- rowMeans(GetAssayData(so, assay="RNA", slot="data")[rownames(so), ])
        res$p_val_adj <- IHW::adj_pvalues(IHW::ihw(res$p_val ~ res$mean_exp, alpha=0.05))

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