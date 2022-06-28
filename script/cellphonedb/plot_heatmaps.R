library(pheatmap)

heatmaps_plot <- function(
    meta_file, 
    pvalues_file, 
    count_network_filename, 
    interaction_count_filename, 
    show_rownames=TRUE, 
    show_colnames=TRUE,    
    scale="none", 
    cluster_cols=FALSE,
    border_color="white", 
    cluster_rows=FALSE, 
    fontsize_row=11,
    fontsize_col=11, 
    treeheight_row=0, 
    family="Arial", 
    treeheight_col=0,
    col1="dodgerblue4", 
    col2="peachpuff", 
    col3="deeppink4", 
    meta_sep="\t", 
    pvalues_sep="\t", 
    pvalue=0.05
){
    ###############
    ### Network ###
    ###############
    
    meta <- read.csv(meta_file, comment.char="", sep=meta_sep)

    all_intr <- read.table(pvalues_file, header=TRUE, stringsAsFactors=FALSE, sep=pvalues_sep, comment.char='', check.names=FALSE)
    intr_pairs <- all_intr$interacting_pair
    all_intr <- all_intr[,-c(1:11)]

    split_sep <- "\\|"
    join_sep <- "|"

    pairs1_all <- unique(meta[,2])

    pairs1 <- c()
    for (i in 1:length(pairs1_all)) {
        
        for (j in 1:length(pairs1_all)) {
          
            pairs1 <- c(pairs1,paste(pairs1_all[i],pairs1_all[j],sep=join_sep))
        
        }
    
    }

    all_count <- matrix(ncol=3)
    colnames(all_count) <- c('SOURCE','TARGET','count')

    count1 <- c()
    for(i in 1:length(pairs1)){
      
        p1 <- strsplit(pairs1[i], split_sep)[[1]][1]
        p2 <- strsplit(pairs1[i], split_sep)[[1]][2]

        n1 <- intr_pairs[which(all_intr[,pairs1[i]]<=pvalue)]
        pairs_rev <- paste(p2, p1, sep=join_sep)
        n2 <- intr_pairs[which(all_intr[ , pairs_rev]<=pvalue)]

        if(p1!=p2) {
          
            count1 <- length(unique(n1))+length(unique(n2))
          
        } else {
          
            count1 <- length(unique(n1))
        }
      
        new_count <- c(p1, p2, count1)
        names(new_count) <- c('SOURCE','TARGET','count')
        all_count <- rbind(all_count, new_count)
    
    }

    all_count <- all_count[-1,]
    
    ##########################
    ### Count interactions ###
    ##########################

    count1 <- c()
    for(i in 1:length(pairs1)) {
      
        p1 <- strsplit(pairs1[i], split_sep)[[1]][1]
        p2 <- strsplit(pairs1[i], split_sep)[[1]][2]

        n1 <- intr_pairs[which(all_intr[,pairs1[i]]<=pvalue)]

        pairs_rev <- paste(p2, p1, sep=join_sep)
        n2 <- intr_pairs[which(all_intr[,pairs_rev]<=pvalue)]
      
        if(p1!=p2) {
          
            count1 <- c(count1,length(unique(n1))+length(unique(n2)))
          
        } else {
          
            count1 <- c(count1,length(unique(n1)))
          
        }
      
  }
  
    if (any(count1)>0) {
        
        count_matrix <- matrix(count1, nrow=length(unique(meta[,2])), ncol=length(unique(meta[,2])))
        rownames(count_matrix) <- unique(meta[,2])
        colnames(count_matrix) <- unique(meta[,2])

        all_sum <- rowSums(count_matrix)
        all_sum <- cbind(names(all_sum), all_sum)

        col.heatmap <- colorRampPalette(c(col1,col2,col3 ))( 1000 )
        
#         order <- c("Ery (1)", "Ery (2)", "Ery (3)", "Ery (4)", "MDP", "MDP Ly6c-", "DC CD11b+ Ly6c+", "DC CD11b+ Esam+", "DC CD11b- Esam+", "DC CD11b- CD8+", "RPM")
#         count_matrix <- count_matrix[order, order]
        
        
        plot <- pheatmap(
            count_matrix, show_rownames=show_rownames, show_colnames=show_colnames, scale=scale, cluster_cols=cluster_cols,
            border_color=border_color, cluster_rows=cluster_rows, fontsize_row=fontsize_row, fontsize_col=fontsize_col,
            main="Interaction count", treeheight_row=treeheight_row, family=family,color=col.heatmap, treeheight_col=treeheight_col, silent=TRUE
        )
        
        plot_log <- pheatmap(
            log10(count_matrix+1), show_rownames=show_rownames, show_colnames=show_colnames, scale=scale, cluster_cols=cluster_cols,
            border_color=border_color, cluster_rows=cluster_rows, fontsize_row=fontsize_row, fontsize_col=fontsize_col,
            main="Interaction log10(count+1)", treeheight_row=treeheight_row, family=family, color=col.heatmap, treeheight_col=treeheight_col, silent=TRUE
        )
    
    } else {
        
        stop("There are no significant results using p-value of: ", pvalue, call.=FALSE)
    
    }
    
    return(list(plot=plot, plot_log=plot_log))
    
}