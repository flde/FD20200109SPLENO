library(ggplot2)
plot_dot <- function(
    selected_rows=NULL,
    selected_columns=NULL,
    means_path='./means.txt',
    pvalues_path='./pvalues.txt',
    means_separator='\t',
    pvalues_separator='\t'
){
    all_pval <-read.table(pvalues_path, header=TRUE, stringsAsFactors=FALSE, sep=means_separator, comment.char='', check.names=FALSE)
    all_means <- read.table(means_path, header=TRUE, stringsAsFactors=FALSE, sep=pvalues_separator, comment.char='', check.names=FALSE)

    intr_pairs <- all_pval$interacting_pair
    all_pval <- all_pval[,-c(1:11)]
    all_means <- all_means[,-c(1:11)]
    
    if(is.null(selected_rows)){selected_rows=intr_pairs}

    if(is.null(selected_columns)){selected_columns=colnames(all_pval)}

    sel_pval <- all_pval[match(selected_rows, intr_pairs), selected_columns]
    sel_means <- all_means[match(selected_rows, intr_pairs), selected_columns]

    df_names <- expand.grid(selected_rows, selected_columns)
    pval <- unlist(sel_pval)
    pval[pval==0] <- 0.0009
    plot_data <- cbind(df_names, pval)
    pr <- unlist(as.data.frame(sel_means))
    pr[pr==0] <- 1
    plot_data <- cbind(plot_data, log2(pr))
    colnames(plot_data) <- c('pair', 'clusters', 'pvalue', 'mean')
    
    plot_data <- plot_data[plot_data$pvalue <= 0.5, ]
    plot_data <- plot_data[abs(plot_data$mean) > 1, ]
    
    color_break <- seq(-(round(max(abs(plot_data$mean)), digits=1) + 0.1), (round(max(abs(plot_data$mean)), digits=1) + 0.1), 0.01)
    color_ramp <- rev(RColorBrewer::brewer.pal(12, "RdBu"))
    my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=399)

    plot <- ggplot(plot_data, aes(x=pair,y=clusters)) +
        geom_point(aes(size=-log10(pvalue), color=mean)) +
        scale_color_gradient2('Log2 mean (Molecule 1, Molecule 2)', midpoint=0, low = "red", mid = "white", high = "blue") +
        theme_bw() +
        theme(
            panel.grid.minor=element_blank(),
            panel.grid.major=element_blank(),
            axis.text=element_text(size=14, colour="black"),
            axis.text.x=element_text(angle=90, hjust=1),
            axis.text.y=element_text(size=12, colour="black"),
            axis.title=element_blank(),
            panel.border=element_rect(size=0.7, linetype="solid", colour="black")
        )

    return(plot)

}
