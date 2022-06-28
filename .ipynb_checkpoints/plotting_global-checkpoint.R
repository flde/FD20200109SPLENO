###############################
### Project global plotting ###
###############################

### Global plotting theme
theme_global_set <- function() {
  
  library(ggplot2, quietly=TRUE)
  
  theme(
    
    panel.grid=element_blank(), 
    panel.background=element_rect(fill=NA, color="black"), 
    
    plot.title=element_text(size=16, face="bold", margin=margin(t=0, r=0, b=5, l=0)), 
#     plot.title=element_blank(), 
    
    axis.title.y=element_text(size=14, face="bold", margin=margin(t=5, r=5, b=5, l=5), angle=90, vjust=2.5), 
    axis.title.x=element_text(size=14, face="bold", margin=margin(t=5, r=5, b=5, l=5)),
    axis.text=element_text(size=12),
    
    legend.text=element_text(size=12), 
    legend.title=element_text(size=12, face="bold"), 
    legend.key=element_blank(), 
    
    strip.text=element_text(size=14, face="bold", margin=margin(t=2.5, r=2.5, b=2.5, l=2.5)), 
    strip.background=element_blank(),
  
  )
  
  }

### Color settings 
library(RColorBrewer, quietly=TRUE)

# CellRanger class
cellranger_class <- c("light grey", RColorBrewer::brewer.pal(9, "Blues")[c(9)])
names(cellranger_class) <- c("Background", "Cell")

# Tissue
tissue <- c(RColorBrewer::brewer.pal(9, "BrBG")[2], RColorBrewer::brewer.pal(9, "BrBG")[8])
names(tissue) <- c("Myeloid", "Progenitor")

# Treatment 
treatment <- c(RColorBrewer::brewer.pal(9, "RdYlBu")[8], RColorBrewer::brewer.pal(9, "RdYlBu")[2])
names(treatment) <- c("NaCl", "CpG")

# Haemosphere main lables 
label_main_haemosphere <- c(
    
    RColorBrewer::brewer.pal(8, "YlGn")[c(4, 8)], # 
    RColorBrewer::brewer.pal(8, "YlOrRd")[1:3],
    RColorBrewer::brewer.pal(8, "Blues")[6:8],
    RColorBrewer::brewer.pal(8, "BuPu")[c(8, 7, 5, 4)],
    RColorBrewer::brewer.pal(8, "RdPu")[6], 
    RColorBrewer::brewer.pal(8, "Reds")[8]

)

names(label_main_haemosphere) <- c(
    
    "MPP", "RPP", 
    "T-cell", "NK", "B-cell", 
    "DC", "Mo", "Mac", 
    "Baso", "Neu", "Eo", "Mast", 
    "Meg", 
    "Ery"
    
)

# Haemosphere fine lables 
label_fine_haemosphere <- c(
    
    RColorBrewer::brewer.pal(8, "YlGn")[3:5], 
    RColorBrewer::brewer.pal(8, "YlGn")[2], 
    RColorBrewer::brewer.pal(8, "YlOrRd")[1:3],
    RColorBrewer::brewer.pal(8, "YlGn")[6], 
    RColorBrewer::brewer.pal(8, "Blues")[4], 
    RColorBrewer::brewer.pal(8, "Blues")[5:8], 
    RColorBrewer::brewer.pal(8, "BuPu")[4:8], 
    RColorBrewer::brewer.pal(8, "RdPu")[4], 
    RColorBrewer::brewer.pal(8, "RdPu")[6], 
    RColorBrewer::brewer.pal(8, "Reds")[4:8], 
    RColorBrewer::brewer.pal(8, "Greys")[5]
    
)

names(label_fine_haemosphere) <- c(
    
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
    "preCFUE", "CFUE", "pbEry", "poEry", "Retic", 
    ## NA 
    "NA"
)

# Cell cycle phase 
cc_phase_class <- c(
    
    RColorBrewer::brewer.pal(8, "Accent")[1:3]

)

names(cc_phase_class) <- c(

    "G1", "S", "G2M"
    
)

leiden_annotation <- c(
    
    RColorBrewer::brewer.pal(9, "PRGn")[1], # Meg
    RColorBrewer::brewer.pal(9, "RdPu")[5], # MEP
    RColorBrewer::brewer.pal(9, "Reds")[6:9], # Ery
    RColorBrewer::brewer.pal(9, "BrBG")[6], # RPM
    RColorBrewer::brewer.pal(9, "BrBG")[7:9], # MDP/Mo
    RColorBrewer::brewer.pal(9, "BrBG")[1:4], # DC
    RColorBrewer::brewer.pal(9, "Blues")[7], # MPP
    RColorBrewer::brewer.pal(9, "PRGn")[7:8], # GMP/Baso 
    RColorBrewer::brewer.pal(9, "PRGn")[3:4] # GMP/Baso
)

names(leiden_annotation) <- c(
    
    "Meg", 
    "MEP",
    "Ery (1)", 
    "Ery (2)", 
    "Ery (3)", 
    "Ery (4)", 
    "RPM", 
    "MDP", 
    "MDP Ly6c-", 
    "Mo",
    "DC CD11b+ Esam+", 
    "DC CD11b- CD8+", 
    "DC CD11b+ Ly6c+",
    "DC CD11b- Esam+",
    "MPP", 
    "GMP",
    "Baso",
    "Plasma cell", 
    "T lymphocyte"

) 
         

# Color list for R
color <- list(
    
    cellranger_class=as.list(cellranger_class), 
    tissue=as.list(tissue), 
    treatment=as.list(treatment), 
    label_main_haemosphere=as.list(label_main_haemosphere), 
    label_fine_haemosphere=as.list(label_fine_haemosphere), 
    leiden_annotation=as.list(leiden_annotation), 
    cc_phase_class=as.list(cc_phase_class)

)


# Get dict for python 
color_dict <- list(
    
    cellranger_class=as.list(cellranger_class), 
    tissue=as.list(tissue), 
    treatment=as.list(treatment), 
    label_main_haemosphere=as.list(label_main_haemosphere), 
    label_fine_haemosphere=as.list(label_fine_haemosphere), 
    leiden_annotation=as.list(leiden_annotation), 
    cc_phase_class=as.list(cc_phase_class)

)

