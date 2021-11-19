###############################
### Project global plotting ###
###############################

### Global plotting theme
theme_global_set <- function() {
  
  library(ggplot2, quietly = TRUE)
  
  theme(
    
    panel.grid = element_blank(), 
    panel.background = element_rect(fill = NA, color = "black"), 
    
    plot.title = element_text(size = 16, face = "bold", margin = margin(t = 0, r = 0, b = 5, l = 0)), 
#     plot.title = element_blank(), 
    
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(t = 2.5, r = 2.5, b = 2.5, l = 2.5), angle = 90, vjust = 2.5), 
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 2.5, r = 2.5, b = 2.5, l = 2.5)),
    axis.text = element_text(size = 10),
    
    legend.text = element_text(size = 10), 
    legend.title = element_text(size = 10, face = "bold"), 
    legend.key = element_blank(), 
    
    strip.text = element_text(size = 12, face = "bold", margin = margin(t = 2.5, r = 2.5, b = 2.5, l = 2.5)), 
    strip.background = element_blank(),
  
  )
  
  }

### Color settings 
library(RColorBrewer, quietly = TRUE)

# CellRanger class
cellranger_class <- c("light grey", RColorBrewer::brewer.pal(9, "Blues")[c(9)])
names(cellranger_class) <- c("Background", "Cell")

# Tissue
tissue <- c(RColorBrewer::brewer.pal(9, "BrBG")[2], RColorBrewer::brewer.pal(9, "BrBG")[8])
names(tissue) <- c("Myeloid", "Progenitor")

# Treatment 
treatment <- c(RColorBrewer::brewer.pal(9, "RdYlBu")[2], RColorBrewer::brewer.pal(9, "RdYlBu")[8])
names(treatment) <- c("CpG", "NaCl")

# Haemosphere main lables 
main_labels_haemosphere <- c(
    
    RColorBrewer::brewer.pal(8, "YlGn")[c(4, 8)], # 
    RColorBrewer::brewer.pal(8, "YlOrRd")[1:3],
    RColorBrewer::brewer.pal(8, "Blues")[6:8],
    RColorBrewer::brewer.pal(8, "BuPu")[c(8, 7, 5, 4)],
    RColorBrewer::brewer.pal(8, "RdPu")[6], 
    RColorBrewer::brewer.pal(8, "Reds")[8]

)

names(main_labels_haemosphere) <- c(
    
    "MPP", "RPP", 
    "T-cell", "NK", "B-cell", 
    "DC", "Mo", "Mac", 
    "Baso", "Neu", "Eo", "Mast", 
    "Meg", 
    "Ery"
    
)

fine_labels_haemosphere <- c(
    
    RColorBrewer::brewer.pal(8, "YlGn")[3:5], 
    RColorBrewer::brewer.pal(8, "YlGn")[2], 
    RColorBrewer::brewer.pal(8, "YlOrRd")[1:3],
    RColorBrewer::brewer.pal(8, "YlGn")[6], 
    RColorBrewer::brewer.pal(8, "Blues")[4], 
    RColorBrewer::brewer.pal(8, "Blues")[5:8], 
    RColorBrewer::brewer.pal(8, "BuPu")[4:8], 
    RColorBrewer::brewer.pal(8, "RdPu")[4], 
    RColorBrewer::brewer.pal(8, "RdPu")[6], 
    RColorBrewer::brewer.pal(8, "Reds")[4:8]
    
)

names(fine_labels_haemosphere) <- c(
    # Multipotent lineage 
    "STHSC", "LSK", "MPP", 
    # Lymphoid lineage 
    "CLP", 
    "T-cell", "NK", "B-cell", 
    # Myeloid lineage 
    "CMP", 
    ## Granulocyte/Monocyte lineage 
    "GMP", 
    "pDC", "cDC", "Mo", "Mac", 
    "Baso", "Neu", "EoP", "Eo", "Mast", 
    ## Megakaryocyte/Erythrocyte lineage 
    "MEP", 
    "Meg", 
    "preCFUE", "CFUE", "pbEry", "poEry", "Retic"
)

cc_phase_class <- c(
    
    RColorBrewer::brewer.pal(8, "Accent")[1:5]

)

names(cc_phase_class) <- c(

    "G1", "G1S", "S", "G2", "G2M"
)


module_class <- c(

    RColorBrewer::brewer.pal(9, "Set1")[1:4], 
    "light gray"

)

names(module_class) <- c(
    
    "Hemoglobin", "Cell cycle", "Mitochondrial", "Ribosomal", "Other"

)

# Get all colors
color <- list(
    
    cellranger_class = cellranger_class, 
    tissue = tissue, 
    treatment = treatment, 
    main_labels_haemosphere = main_labels_haemosphere, 
    fine_labels_haemosphere = fine_labels_haemosphere, 
    cc_phase_class = cc_phase_class, 
    module_class = module_class

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
