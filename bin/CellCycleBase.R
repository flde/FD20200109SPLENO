####################
### cc_whitfield ###
####################

cc_genes_whitfield_load <- function() {
    
    # Set ssl for biomart  
    httr::set_config(httr::config(ssl_verifypeer=FALSE))
    
    cc_whitfield <- read.csv("data/cell_cycle/cc_whitfield.csv")

    # Select and rename columns
    cc_whitfield <- unique(cc_whitfield[, c("LLID", "PHASE")])
    colnames(cc_whitfield) <- c("entrezgene_id", "cc_phase")
    # Use only unique entrenzgene_id
    cc_whitfield <- cc_whitfield[!cc_whitfield$entrezgene_id %in% cc_whitfield$entrezgene_id[duplicated(cc_whitfield$entrezgene_id)], ]
    # Select cell cycle phase of interest 
    cc_whitfield$cc_phase <- gsub("S phase", "S", cc_whitfield$cc_phase)
    cc_whitfield <- cc_whitfield[cc_whitfield$cc_phase %in% c("G1/S", "S", "G2", "G2/M"), ]

    # Convert human entrenz gene ID to human symbol
    library(biomaRt) 
    ensembl_hgnc <- useMart("ensembl", dataset="hsapiens_gene_ensembl") 
    symbol_hgnc <- getBM(attributes=c("entrezgene_id", "hgnc_symbol"), filters="entrezgene_id", values=cc_whitfield$entrezgene_id, mart=ensembl_hgnc)

    # Convert human symbol to mouse symbol 
    ensembl_mm <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
    symbol_mm <- getLDS(attributes=c("hgnc_symbol"), filters="hgnc_symbol", values=symbol_hgnc$hgnc_symbol, mart=ensembl_hgnc, attributesL=c("mgi_symbol"), martL=ensembl_mm, uniqueRows=T)
    colnames(symbol_mm) <- c("hgnc_symbol", "mgi_symbol")

    # Combine and filter conversion results 
    convert <- dplyr::left_join(symbol_hgnc, symbol_mm, by="hgnc_symbol")
    convert <- convert[!convert$entrezgene_id %in% convert$entrezgene_id[duplicated(convert$entrezgene_id)], ]
    convert <- convert[!convert$mgi_symbol %in% convert$mgi_symbol[duplicated(convert$mgi_symbol)], ]

    # Add mgi_symbol to whitfield cell cycle phase annotation 
    cc_whitfield <- dplyr::left_join(cc_whitfield, convert, by="entrezgene_id")
    # Remove entries without mouse gene symbol annotation 
    cc_whitfield <- na.omit(cc_whitfield)
    # cc_whitfield <- cc_whitfield[cc_whitfield$mgi_symbol %in% rownames(so), ]

    # Combine to list 
    cc_genes_whitfield <- list(
        g1s=cc_whitfield[cc_whitfield$cc_phase == "G1/S", ]$mgi_symbol,
        s=cc_whitfield[cc_whitfield$cc_phase == "S", ]$mgi_symbol,
        g2=cc_whitfield[cc_whitfield$cc_phase == "G2", ]$mgi_symbol,
        g2m=cc_whitfield[cc_whitfield$cc_phase == "G2/M", ]$mgi_symbol
    )
    
    # Transform cc_genes list to data frame 
    cc_genes_whitfield <- melt(cc_genes_whitfield)
    colnames(cc_genes_whitfield) <- c("gene", "cc_phase_class")

    return(cc_genes_whitfield)
    
}

#################
### cc_seurat ###
#################

cc_genes_seurat_load <- function() {
    
    # Set ssl for biomart  
    httr::set_config(httr::config(ssl_verifypeer=FALSE))
    
    # Get human cell cycle genes from seurat 
    require(Seurat)
    cc_genes_seurat_s <- cc.genes.updated.2019$s.genes
    cc_genes_seurat_g2m <- cc.genes.updated.2019$g2m.genes
    
    # Load Biomart 
    library(biomaRt) 
    ensembl_hgnc <- useMart("ensembl", dataset="hsapiens_gene_ensembl") 
    ensembl_mm <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

    # Translate human to mouse 
    cc_genes_seurat_s <- getLDS(attributes=c("hgnc_symbol"), filters="hgnc_symbol", values=cc_genes_seurat_s, mart=ensembl_hgnc, attributesL=c("mgi_symbol"), martL=ensembl_mm, uniqueRows=T)
    cc_genes_seurat_s <- cc_genes_seurat_s[, 2]

    cc_genes_seurat_g2m <- getLDS(attributes=c("hgnc_symbol"), filters="hgnc_symbol", values=cc_genes_seurat_g2m, mart=ensembl_hgnc, attributesL=c("mgi_symbol"), martL=ensembl_mm, uniqueRows=T)
    cc_genes_seurat_g2m <- cc_genes_seurat_g2m[, 2]

    # Combine to list 
    cc_genes_seurat <- list(s=cc_genes_seurat_s, g2m=cc_genes_seurat_g2m)
    
    # Transform cc_genes to data frame 
    cc_genes_seurat <- melt(cc_genes_seurat)
    colnames(cc_genes_seurat) <- c("gene", "cc_phase_class")
    
    return(cc_genes_seurat)
    
}

