library_load <- suppressMessages(
    
    suppressWarnings(
        
        list(

            # Reactome GSEA/ORA and data base 
            library(ReactomePA),
            library(clusterProfiler), # Symbol to ID
            library(org.Mm.eg.db),

            library(AnnotationDbi),
            library(reactome.db),

            # Data 
            library(dplyr), 
            library(tidyverse), 
            library(data.table)

        )
    )
)

#################################
### Get ranks from DEA result ###
#################################
dea_ranks <- function(dea_res) {

    # Make ranks 
    dea_res$p_val_adj <- ifelse(dea_res$p_val_adj==0, min(dea_res$p_val_adj[dea_res$p_val_adj>0]), dea_res$p_val_adj)
    dea_res$sign_log_adj_p_values <- -log10(dea_res$p_val_adj) * sign(dea_res$avg_log2FC)

    ranks <- dea_res$sign_log_adj_p_values
    names(ranks) <- dea_res$gene
    ranks <- ranks[order(ranks)]
    ranks <- rev(ranks)

    

    ranks <- ranks[abs(ranks)>0]

    return(ranks)
    
}

########################
### GSEA Reactome PA ### 
########################
gsea_pa <- function(ranks, convert_names=FALSE) {

    # Convert gene symbol to ID 
    if(convert_names) {
        
        names(ranks) <- clusterProfiler::bitr(names(ranks), fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Mm.eg.db, drop=FALSE)$ENTREZID
        ranks <- ranks[!is.na(names(ranks))]
        
    }

    pt_gsea_res <- ReactomePA::gsePathway(ranks, organism="mouse", pvalueCutoff=1, minGSSize=5, maxGSSize=500, pAdjustMethod="BH") %>% as.data.frame()

    return(pt_gsea_res)
    
}

#################################
### Reactome hierarchy parser ###
#################################
parse_hierarchy <- function(enrichment_res, level=3) {

    # Load Reactome hierarchy
    rh <- readRDS("data/reference/reactome/reactome_hierarchy.rds")

    if(level==3) {

        # Get IDs from higher level that are not terminal levels 
        enrichment_res <- enrichment_res %>% dplyr::filter(!(ID %in% rh[!is.na(rh$ID_2), ]$ID_1 | ID %in% rh[!is.na(rh$ID_3), ]$ID_2))
    
        # Filter for level three except there is no such level in the hierarchy 
        enrichment_res <- enrichment_res %>% dplyr::filter((ID %in% rh$ID_3) | (ID %in% rh[is.na(rh$ID_2), ]$ID_1 | ID %in% rh[is.na(rh$ID_3), ]$ID_2))  
        
    }


    if(level==4) {

        # Get IDs from higher level that are not terminal levels 
        enrichment_res <- enrichment_res %>% dplyr::filter(!(ID %in% rh[!is.na(rh$ID_2), ]$ID_1 | ID %in% rh[!is.na(rh$ID_3), ]$ID_2 | ID %in% rh[!is.na(rh$ID_4), ]$ID_3))
    
        # Filter for level four except there is no such level in the hierarchy 
        enrichment_res <- enrichment_res %>% dplyr::filter((ID %in% rh$ID_4) | (ID %in% rh[is.na(rh$ID_2), ]$ID_1 | ID %in% rh[is.na(rh$ID_3), ]$ID_2 | ID %in% rh[is.na(rh$ID_4), ]$ID_3))  
        
    }

    # Set BH for levels 
    enrichment_res$p.adjust=p.adjust(enrichment_res$pvalue, method="BH")

    # Annotate level 1
    enrichment_res <- lapply(enrichment_res$ID, function(i) {data.frame(ID=i, top_level=rh[apply(rh, 1, function(x) any(grepl(i, x))), ]$pathway_1) %>% dplyr::distinct()}) %>% do.call("rbind", .) %>% dplyr::left_join(enrichment_res, ., by=join_by(ID)) %>% dplyr::distinct()

    return(enrichment_res)
    
}
######################################
### Get Reactome gene set by level ###
######################################
reactome_gs <- function(rh_level=3) {

    rh <- readRDS("data/reference/reactome/reactome_hierarchy.rds")

    reactome <- rh[, c(paste0("ID_", rh_level), paste0("pathway_", rh_level))] %>% na.omit() %>% dplyr::distinct()
    colnames(reactome) <- c("PATHID", "PATHNAME")
    
    map_1 <- AnnotationDbi::select(
                
        reactome.db,
        keys=reactome$PATHID,
        keytype="PATHID",
        columns="ENTREZID"
        
    )
    
    map_2 <- AnnotationDbi::select(
                
        org.Mm.eg.db,
        keys=map_1$ENTREZID,
        keytype="ENTREZID",
        columns="SYMBOL"
        )
    
    gene_set <- dplyr::left_join(map_1, map_2, by=join_by(ENTREZID)) %>% dplyr::left_join(reactome, ., by=join_by(PATHID)) %>% dplyr::select(-ENTREZID, -PATHID) %>% na.omit() %>% dplyr::distinct() %>% dplyr::rename(gs_name=PATHNAME, gene_symbol=SYMBOL)

    return(gene_set)
}
                                                                                           
################################
### Reactome top level order ### 
################################
top_level_order <- c("Cell Cycle", "DNA Repair", "DNA Replication", "Metabolism", "Metabolism of RNA", "Metabolism of proteins", "Gene expression (Transcription)", "Cellular responses to stimuli", "Signal Transduction", "Immune System", "Programmed Cell Death", "Autophagy", "Vesicle-mediated transport", "Drug ADME", "Other")

####################
### PT GSEA plot ###
####################
pt_gsea_pl_pa <- function(pt_gsea_res, adj_pval_thr=0.1, title=NULL, color=c(RColorBrewer::brewer.pal(8, "Set1")[1], RColorBrewer::brewer.pal(8, "Set1")[2]), color_names=c("Pos", "Neg"), size_range=5, top=NULL) {

    # Fix names 
    pt_gsea_res <- pt_gsea_res %>% dplyr::rename("pval"="pvalue", "padj"="p.adjust", "pathway"="Description")

    # Exclude terms 
    pt_gsea_res <- pt_gsea_res %>% dplyr::filter(!pathway %in% grep("Neutrophil|B\\scell|T\\scell|TCR|BCR", pt_gsea_res$pathway, value=TRUE))

    # Set GSEA data frame 
    pt_gsea_res <- as.data.frame(pt_gsea_res)
    pt_gsea_res <- na.omit(pt_gsea_res) 
    
    # Set color names 
    color <- setNames(color, color_names)
    
    # Filter hits 
    gsea_up <- pt_gsea_res[sign(pt_gsea_res$NES)==+1, ]
    gsea_down <- pt_gsea_res[sign(pt_gsea_res$NES)==-1, ]

    if(!is.null(top)) {
        
        gsea_up <- gsea_up[order(gsea_up$padj), ][1:top, ]
        gsea_down <- gsea_down[order(gsea_down$padj), ][1:top, ]
        
    } else {

        gsea_up <- gsea_up[gsea_up$padj <= adj_pval_thr, ]
        gsea_down <- gsea_down[gsea_down$padj <= adj_pval_thr, ]

    }

    pt_gsea_res <- rbind(gsea_up, gsea_down)
    pt_gsea_res <- na.omit(pt_gsea_res)
    pt_gsea_res <- distinct(pt_gsea_res)
    
    # Add color 
    pt_gsea_res$color <- ifelse(sign(pt_gsea_res$NES)==1, names(color)[1], names(color)[2])
    pt_gsea_res$color <- ifelse(pt_gsea_res$padj<=adj_pval_thr, pt_gsea_res$color, NA)

    pt_gsea_res$sign_log_pval_values <- -log10(pt_gsea_res$padj) * sign(pt_gsea_res$NES)
    
    # Order hits 
    pt_gsea_res <- pt_gsea_res[rev(order(pt_gsea_res$sign_log_pval_values)), ]

    if("top_level" %in% colnames(pt_gsea_res)) {

        pt_gsea_res <- pt_gsea_res %>% dplyr::group_by(top_level) %>% dplyr::arrange(sign_log_pval_values) %>% dplyr::ungroup() %>% 
            dplyr::mutate(pathway=interaction(top_level, pathway, sep="_")) %>% 
            dplyr::mutate(pathway=factor(pathway, levels=rev(pathway)))
        
    } else {

        pt_gsea_res$pathway <- factor(pt_gsea_res$pathway, levels=rev(pt_gsea_res$pathway))
        
    }
    
    x_max <- max(abs(pt_gsea_res$sign_log_pval_values))
    if(x_max<abs(log10(adj_pval_thr))) {x_max <- abs(log10(adj_pval_thr))}
    x_max <- ceiling(ceiling(x_max))
    
    # Plot 
    plot <- ggplot(pt_gsea_res, aes(x=sign_log_pval_values, y=pathway, color=color, size=abs(NES))) + 
        
        geom_vline(xintercept=-log10(adj_pval_thr), linetype="dashed") + 
        geom_vline(xintercept=log10(adj_pval_thr), linetype="dashed") +
    
        geom_point() +

        ggtitle(title) +
        xlab("Signed -log10 adj. p-value") + ylab("") + 
        scale_x_continuous(breaks=c(-x_max, 0, x_max), limits=c(-x_max, x_max)) +
        scale_y_discrete(labels = function(x) sub(".*_", "", x)) + 
        scale_color_manual(values=color, na.value="black", drop=FALSE) +
        scale_size_continuous(limits=c(0, ceiling(max(abs(pt_gsea_res$NES)))), breaks=seq(0, ceiling(max(abs(pt_gsea_res$NES)))), range=c(1, 5)) +
        guides(
            
            color=guide_legend(order=1, title="Agent", keywidth=0.75, keyheight=0.75, override.aes=list(size=4)), 
            size=guide_legend(order=2, title="Abs. (NES)", keywidth=0.75, keyheight=0.75)
            
        ) +
    
        theme(
            
            legend.position="bottom", 
            legend.justification="top", 
            axis.text.y=element_text(size=14, hjust=1, vjust=0.5, face="plain", margin=margin(t=0, r=2, b=0, l=0), color="black")
            
        ) 
    
    if("top_level" %in% colnames(pt_gsea_res)) {plot <- plot + ggforce::facet_col(~top_level, scales="free_y", space="free")}

    return(plot)
    
}

