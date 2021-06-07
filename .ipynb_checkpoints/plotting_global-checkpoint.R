###############################
### Project global plotting ###
###############################

### Global plotting theme
theme_global_set <- function() {
  
  require(ggplot2)
  
  theme(
    
    panel.grid = element_blank(), 
    panel.background = element_rect(fill = NA, color = "black"), 
    
    # plot.title = element_text(size = 12, face = "bold", margin = margin(t = 0, r = 0, b = 5, l = 0)), 
    plot.title = element_blank(), 
    
    axis.title.y = element_text(size = 10, face = "bold", margin = margin(t = 2.5, r = 2.5, b = 2.5, l = 2.5), angle = 90, vjust = 2.5), 
    axis.title.x = element_text(size = 10, face = "bold", margin = margin(t = 2.5, r = 2.5, b = 2.5, l = 2.5)),
    axis.text = element_text(size = 8),
    
    legend.text = element_text(size = 8), 
    legend.title = element_text(size = 8, face = "bold"), 
    legend.key = element_blank(), 
    
    strip.text = element_text(size = 10, face = "bold", margin = margin(t = 2.5, r = 2.5, b = 2.5, l = 2.5)), 
    strip.background = element_blank(),
  
  )
  
  }

# Color settings 
require(RColorBrewer)

cellranger_class <- c("light grey", RColorBrewer::brewer.pal(9, "Blues")[c(9)])
names(cellranger_class) <- c("Background", "Cell")

tissue <- c(RColorBrewer::brewer.pal(9, "BrBG")[2], RColorBrewer::brewer.pal(9, "BrBG")[8])
names(tissue) <- c("Myeloid", "Progenitor")

treatment <- c(RColorBrewer::brewer.pal(9, "RdYlBu")[2], RColorBrewer::brewer.pal(9, "RdYlBu")[8])
names(treatment) <- c("CpG", "NaCl")

main_labels <- c(
  
  RColorBrewer::brewer.pal(9, "Purples")[6:9], 
  RColorBrewer::brewer.pal(9, "Greens")[5:6], 
  RColorBrewer::brewer.pal(9, "PiYG")[1:3], 
  RColorBrewer::brewer.pal(9, "PiYG")[4],
  RColorBrewer::brewer.pal(9, "PuBu")[6:9], 
  RColorBrewer::brewer.pal(9, "BrBG")[8], 
  RColorBrewer::brewer.pal(9, "BrBG")[8:9]
  
  )

names(main_labels) <- c(
  
  "T cells", "NK cells", "NKT", "Tgd",
  "B cells", "B cells, pro", 
  "Basophils", "Neutrophils", "Eosinophils", 
  "Mast cells", 
  "Monocytes", "Macrophages", "DC", "ILC",
  "Stem cells", 
  "Endothelial cells", "Epithelial cells"
  
)

so_color <- list(
    
    cellranger_class = cellranger_class, 
    tissue = tissue, 
    treatment = treatment, 
    main_labels = main_labels

)

# Integer breaks 
integer_breaks <- function(n = 5, ...) {
  fxn <- function(x) {
    breaks <- floor(pretty(x, n, ...))
    names(breaks) <- attr(breaks, "labels")
    breaks
  }
  return(fxn)
}
