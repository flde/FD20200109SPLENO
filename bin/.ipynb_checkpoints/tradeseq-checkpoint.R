#########################
### TradeSeq workflow ###
#########################
tradeseq_workflow <- function(so, pseudotime, cell_weights=NULL, conditions=NULL, genes=NULL, suffix="", log2_thr=0, nknots=6, parallel=TRUE, workers=16, compute=FALSE) {
    
    if(compute) {
        
        # Prallel computing 
        print(future::availableCores())
        BPPARAM <- BiocParallel::bpparam()
        BPPARAM <- MulticoreParam(workers=future::availableCores())        
        
        # Set cell weights 
        if(is.null(cell_weights)) {cell_weights <- rep(1, ncol(so))}
        
        # Get and filter matrix 
        counts <- GetAssayData(so, assay="RNA", slot="counts")
        
        # Filter count matrix by genes
        if(!is.null(genes)) {counts <- counts[rownames(counts) %in% genes, ]}
        
        # evaluateK
        evaluate_k <- evaluateK(
            
            counts=counts,
            pseudotime=pseudotime,
            cellWeights=cell_weights,
            conditions=conditions,
            nGenes=500,
            k=3:15, 
            parallel=parallel, 
            BPPARAM=BPPARAM

        )
        
        # fitGAM
        fitgam <- fitGAM(

            counts=counts, 
            pseudotime=pseudotime,
            cellWeights=cell_weights,
            conditions=conditions, 
            genes=genes, 
            nknots=nknots,
            verbose=TRUE,
            parallel=parallel, 
            BPPARAM=BPPARAM

        )
        
        # associationTest
        association <- associationTest(fitgam, l2fc=log2_thr, lineage=TRUE, contrastType="end", inverse="Chol")

        # Save results 
        tradeseq_result <- list(fitgam=fitgam, association=association, evaluate_k=evaluate_k)
        
        saveRDS(tradeseq_result, paste0("result/tradeseq/tradeseq", suffix, ".rds"))
        
        
    } else {
        
        tradeseq_result <- readRDS(paste0("result/tradeseq/tradeseq", suffix, ".rds"))
        
    }
    
    return(tradeseq_result)

}

################
### Heat map ###
################
hm_smooth <- function(fitgam, genes, condition_test=NULL, n, width=1.0, height=0.005, padj_thr=0.05, cluster_rows=FALSE, show_row_names=TRUE, user_raster=FALSE) {
    
    mat <- predictSmooth(fitgam, gene=genes, nPoints=n, tidy=FALSE) 
    mat <- t(apply(mat, 1, scales::rescale))
    
    breaks <- seq(0, 1, length.out=4)
    color_ramp <- c("white", rev(brewer.pal(11,"RdBu"))[9:11])
    color_ramp <- viridis::mako(length(breaks))

    top_annotation <- HeatmapAnnotation(

        df=data.frame(treatment=factor(c(rep("NaCl", n), rep("CpG", n)), levels=c("NaCl", "CpG"))), 
        col=list(treatment=unlist(color$treatment)), 
        simple_anno_size=unit(4, "mm"), 
        show_legend=FALSE

    )

    row_order <- pheatmap(

        mat,
        cluster_cols=FALSE,
        cluster_rows=TRUE, 
        silent=TRUE

    )
    
    if(!is.null(condition_test)) {
        
        # Left row annotation 
        condition_test$padj <- p.adjust(condition_test$pvalue, "fdr")
        condition_test$padj <- ifelse(condition_test$padj==0, min(condition_test$padj[condition_test$padj > 0], na.rm=TRUE), condition_test$padj)
        
        condition_test <- condition_test[genes, ]

        padj_annotation <- condition_test$padj
        padj_annotation[padj_annotation>padj_thr] <- NA
        padj_annotation <- -log10(padj_annotation)

        padj_annotation_color <- brewer.pal(9,"BuPu")[c(1, 9)]
        padj_annotation_breaks <- c(0, max(padj_annotation, na.rm=TRUE))

        if(is.infinite(padj_annotation_breaks[2])) {padj_annotation_breaks[2] <- -log10(0.1)} # If all padj NA the max value gets Inf. I set a high pval instead 
        
        row_annotation <- HeatmapAnnotation(
            
            df=data.frame(padj_annotation=padj_annotation), 
            annotation_legend_param=list(padj_annotation=list(title="-log10(padj)", title_gp=gpar(fontsize=20, fontface="plain"), labels_gp=gpar(fontsize=18))), 
            col=list(padj_annotation=colorRamp2(padj_annotation_breaks, padj_annotation_color)), 
            gp=gpar(col="white", lwd=1), 
            na_col="white", 
            which="row", 
            height=unit(5.5, "in"),
            width=unit(0.75, "in"), 
            show_annotation_name=FALSE
             
         )
        
    } else {
        
        row_annotation=NULL
        
    }
   
    hm_1 <- Heatmap(

        matrix=mat,   
        name="z-score", 
        column_title_gp=gpar(fontsize=20, fontface="plain"), 
        column_names_gp =grid::gpar(fontsize=20, fontface="plain"), 
        row_title_gp=gpar(fontsize=18, fontface="plain"),
        row_names_gp=grid::gpar(fontsize=18, fontface="plain"), 
        col=colorRamp2(breaks, color_ramp),
        cluster_rows=cluster_rows, 
        cluster_columns=FALSE,
        show_row_names=show_row_names,
        show_column_names=FALSE, 
        top_annotation=top_annotation, 
        left_annotation=row_annotation, 
        row_dend_width=unit(0, "cm"), 
        width=n*unit(width, "mm"), 
        height=nrow(mat)*unit(height, "mm"),
        show_row_dend=FALSE, 
        column_split=factor(c(rep("NaCl", n), rep("CpG", n)), levels=c("NaCl", "CpG")), 
        cluster_row_slices=FALSE, 
        heatmap_legend_param=list(title_gp=gpar(fontsize=20, fontface="plain"), labels_gp=gpar(fontsize=18), legend_height=unit(10, "mm"), grid_width=unit(4, "mm")), 
        use_raster=user_raster, 
        raster_resize_mat=mean

    ) %>% as.ggplot()

    return(hm_1)  
    
    
}
#######################
### Import cellrank ###
#######################
import_cellrank <- function(so, suffix="_eb", absorption_probability_col=NULL) {
    
    # DPT 
    dpt_pseudotime <- read.csv(paste0("result/cellrank/dpt_pseudotime", suffix, ".csv"), row.names=1)
    so <- subset(so, subset=cell_id %in% rownames(dpt_pseudotime))
    so <- AddMetaData(so, dpt_pseudotime)
    
    # DPT treatment 
    dpt_pseudotime_nacl <- read.csv(paste0("result/cellrank/dpt_pseudotime", suffix, "_nacl", ".csv"), row.names=1)
    dpt_pseudotime_cpg <- read.csv(paste0("result/cellrank/dpt_pseudotime", suffix, "_cpg", ".csv"), row.names=1)
    dpt_pseudotime_treatment <- rbind(dpt_pseudotime_nacl, dpt_pseudotime_cpg)
    colnames(dpt_pseudotime_treatment) <- "dpt_pseudotime_treatment"
    so <- AddMetaData(so, dpt_pseudotime_treatment)
    
    # FA graph 
    so[["fag"]] <- CreateDimReducObject(embeddings=as.matrix(read.csv(paste0("result/cellrank/fag", suffix, ".csv"), row.names=1)), key="UMAP_")
    
    # convert back to singleCellExperiment
    sce <- as.SingleCellExperiment(so, assay="RNA")
        
    # Condiments imbalance score  
    imbalance_score_eb <- condiments::imbalance_score(Object=so@reductions$umap@cell.embeddings, conditions=so$treatment, k=20, smooth=40)
    so$imbalance_score <- imbalance_score_eb$scaled_scores
    
    # Plotting 
    dplot_1 <- dplot(so, reduction="fag", group_by="cell_type_fine_detail", alpha=1) + scale_color_manual(values=color[["cell_type_fine_detail"]][names(color[["cell_type_fine_detail"]]) %in% so$cell_type_fine]) 
    dplot_2 <- dplot(so, reduction="fag", group_by="treatment", alpha=1) + scale_color_manual(values=color[["treatment"]])
    fplot_1 <- fplot(so, reduction="fag", features="EB..5.") + ggtitle("Absorption probability") + scale_color_viridis(option="G")
    fplot_2 <- fplot(so, reduction="fag", features="imbalance_score") + ggtitle("Imbalance score") + scale_color_viridis(option="G")
    fplot_3 <- fplot(so, reduction="fag", features="dpt_pseudotime") + ggtitle("DPT pseudotime") + scale_color_viridis(option="G") 
    fplot_4 <- fplot(so, reduction="fag", features="dpt_pseudotime_treatment") + ggtitle("DPT pseudotime treatment") + scale_color_viridis(option="G") 

    print(dplot_1 + dplot_2 + fplot_1 + fplot_2 + fplot_3 + fplot_4 + plot_layout(ncol=6, nrow=1))
    
    return(so)
    
}

##########################
### Pseudotime density ###
##########################
pseudotime_density_plot <- function(so, pseudotime="dpt_pseudotime"){
    
    plot <- ggplot(so@meta.data, aes_string(x=pseudotime)) + 
        geom_density(alpha=0.8, aes(fill=treatment), color="transparent") + 
        geom_density(aes(col=treatment), fill="transparent", guide=FALSE, size=1.5) +
        labs(x="Pseudotime", fill="Treatment") +
        guides(color="none") +
        scale_color_manual(values=color[["treatment"]]) +
        scale_fill_manual(values=color[["treatment"]]) 

    return(plot)

}

######################
### Condition GSEA ###
######################
gsea_condition <- function(stats, species="Mus musculus", category="C5", subcategory="BP") {
    
    # Get gene set
    gene_set <- msigdbr(species=species, category=category, subcategory=subcategory) 
    gene_set <- split(gene_set, x=gene_set$gene_symbol, f=gene_set$gs_name)

    # Retain only pathways with at least three genes from condition test
    gene_set_filter <- lapply(gene_set, function(x) {sum(names(stats) %in% x)>=3})
    gene_set <- gene_set[unlist(gene_set_filter)]

    # Filter gene sets by number of genes
    gene_set_filter <- lapply(gene_set, function(x) {length(x)>=5})
    gene_set <- gene_set[unlist(gene_set_filter)]
    
    # Run fast GSEA
    gsea_result <- fgsea(pathways=gene_set, stats=stats, nperm=10^5, minSize=5)
    
    return(gsea_result)
    
}

###########################
### Condition GSEA plot ###
###########################
gsea_plot <- function(gsea, pval_thr=0.05, title="", pathway_suffix=NULL, pathway_filter=NULL, scale_size=2, top=10) {
    
    gsea <- as.data.frame(gsea)
    
    gsea$sig_class <- ifelse(gsea$pval<=pval_thr, paste0("pval<=", pval_thr), NA)

    gsea$minus_log_pval <- -log10(gsea$pval)
    gsea <- arrange(gsea, pval)[1:top, ]

    if(!is.null(pathway_filter)) {gsea <- gsea[gsea$pathway %in% pathway_filter, ]}
    if(!is.null(pathway_suffix)) {gsea$pathway <- gsub(pathway_suffix, "", gsea$pathway)}
    gsea$pathway <- gsub("_", " ", gsea$pathway)
    
    # Order hits 
    gsea$pathway <- factor(gsea$pathway, levels=rev(gsea$pathway))

    x_max <- max(gsea$minus_log_pval)
    
    plot <- ggplot(gsea, aes(x=minus_log_pval, y=pathway, title=title, size=abs(NES), color=sig_class)) + 
        geom_point() + 
        ggtitle(title) + xlab("-log10(pval.)") + ylab("") +
        scale_x_continuous(expand=c(0.1, 0.1)) + 
        scale_y_discrete(expand=c(0, 1)) + 
        scale_color_manual(values=brewer.pal(11,"RdBu")[2], na.value="grey") + 
        scale_size(range=c(0, scale_size)) + 
        guides(
            
            color=guide_legend(title="Significance", keywidth=0.75, keyheight=0.75), 
            size=guide_legend(order=2, title="Abs. (NES)", keywidth=0.75, keyheight=0.75)
        
        ) + 
    
        theme(axis.text.y=element_text(hjust=1, vjust=0.5))
    
    return(plot)
    
}

###############################
### Plort smooth expression ###
###############################
plot_smooth <- function(sce, genes, point=TRUE, line=TRUE) {
    
    cnt <- assays(sce)$counts[genes, , drop=FALSE]
    cnt <- as.data.frame(t(cnt))
    colnames(cnt) <- "exp"

    cnt$pseudotime <- colData(sce)$crv$pseudotime
    cnt$condition <- colData(sce)$tradeSeq$conditions
    
    cnt_smooth <- predictSmooth(sce, genes, nPoints=50, tidy=TRUE)
    cnt_smooth <- cnt_smooth[, c("yhat", "time", "condition")]
    colnames(cnt_smooth) <- c("exp", "pseudotime", "condition")
    
    p <- ggplot(NULL, aes(x=pseudotime, y=log1p(exp), color=condition)) + 
        scale_color_manual(values=unlist(color$treatment))
    
    if(point) {p <- p + geom_point(data=cnt, size=1)}
    if(line) {p <- p + geom_line(data=cnt_smooth, size=3, alpha=1)}
        
    
    return(p)
    
}
