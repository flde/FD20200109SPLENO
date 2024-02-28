###############################
### Project global plotting ###
###############################
theme_global_set <- function(size_select=2) {
    
    size <- data.frame(
    
        size_1=c(16, 18, 20),
        size_2=c(10, 12, 16),
        size_3=c(8, 10, 12)
        
    )

    size <- size[, size_select]
    
    library(ggplot2, quietly=TRUE)
  
    theme(
        
        panel.background=element_blank(), 
        panel.spacing=unit(0.2, "lines"),

        plot.title=element_text(size=size[3], face="plain", margin=margin(t=0, r=0, b=10, l=0), color="black"), 
    
        axis.title.y=element_text(size=size[2], face="plain", margin=margin(t=0, r=5, b=0, l=0), angle=90, color="black"), 
        axis.title.x=element_text(size=size[2], face="plain", margin=margin(t=5, r=0, b=0, l=0), color="black"),
        axis.text.y=element_text(size=size[1], face="plain", margin=margin(t=0, r=2, b=0, l=0), color="black"),
        axis.text.x=element_text(size=size[1], face="plain", margin=margin(t=2, r=0, b=0, l=0), color="black"),
      
        axis.ticks=element_line(color="black", size=unit(1/2.141959, "pt")), 
        axis.line=element_line(color="black", size=unit(1/2.141959, "pt")), 
      
        legend.key.height=unit(0.5, "cm"), 
        legend.key.width=unit(0.5, "cm"),
        legend.key.size=unit(0.5, 'cm'),
      
        legend.title=element_text(size=size[2], face="plain", color="black"),
        legend.text=element_text(size=size[1], face="plain", color="black"), 
        legend.key=element_rect(fill = "transparent", colour = "transparent"), 
    
        strip.text=element_text(size=size[2], margin=margin(t=2, r=2, b=2, l=2), face="plain", color="black") 
    
    )

}

########################################
### Remove layers from ggplot object ###
########################################
remove_geom <- function(ggplot2_object, geom_type) {
    
  # Delete layers that match the requested type.
  layers <- lapply(ggplot2_object$layers, function(x) {
      
    if (class(x$geom)[1] == geom_type) {
        
      NULL
        
    } else {
      x
    }
  })
    
  # Delete the unwanted layers.
  layers <- layers[!sapply(layers, is.null)]
  ggplot2_object$layers <- layers
  return(ggplot2_object)
    
}

######################
### Color settings ###
######################

library(RColorBrewer, quietly=TRUE)

# CellRanger class
cellranger_class <- c("light grey", RColorBrewer::brewer.pal(9, "Blues")[c(9)])
names(cellranger_class) <- c("Background", "Cell")

# Tissue
tissue <- c(RColorBrewer::brewer.pal(9, "BrBG")[2], RColorBrewer::brewer.pal(9, "BrBG")[8])
names(tissue) <- c("Myeloid", "Progenitor")

# Treatment 
treatment <- c("#7F7F7F", "#722B90")
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
    RColorBrewer::brewer.pal(8, "Blues")[3:4], 
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
    "T-cell", 
    "NK", 
    "B-cell", 
    # Myeloid lineage 
    "EoP", 
    ## Granulocyte/Monocyte lineage 
    "GMP", 
    "CMP", 
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

cell_type_main <- c(
    
    RColorBrewer::brewer.pal(9, "YlOrRd")[2], # MLP
    RColorBrewer::brewer.pal(9, "BrBG")[7], # GMP
    RColorBrewer::brewer.pal(9, "Blues")[5:7], # MDP, cMo, ncMo
    rev(RColorBrewer::brewer.pal(9, "BrBG")[c(3,2)]), # RPM
    RColorBrewer::brewer.pal(9, "PRGn")[2], # cDC1
    RColorBrewer::brewer.pal(10, "PRGn")[9], # cDC2
    RColorBrewer::brewer.pal(9, "PiYG")[2], # MegP
    colorRampPalette(c(RColorBrewer::brewer.pal(9, "OrRd")[3], RColorBrewer::brewer.pal(9, "OrRd")[9]))(13)[3], # MEP
    colorRampPalette(c(RColorBrewer::brewer.pal(9, "OrRd")[3], RColorBrewer::brewer.pal(9, "OrRd")[9]))(13)[7], # ProEB
    colorRampPalette(c(RColorBrewer::brewer.pal(9, "OrRd")[3], RColorBrewer::brewer.pal(9, "OrRd")[9]))(13)[12] # EB
    
    
)

names(cell_type_main) <- c(
        
    "MLP",

    "GMP",

    "MDP",
    "cMo", 
    "ncMo", 

    "RPM",  
    "PreRPM", 

    "cDC1", 
    
    "cDC2",
    
    "MegP",
    "MEP",
    "ProEB",
    "EB"

) 

cell_type_fine <- c(
    
    RColorBrewer::brewer.pal(9, "YlOrRd")[2], # MLP
    RColorBrewer::brewer.pal(9, "BrBG")[7:9], # GMP
    RColorBrewer::brewer.pal(9, "Blues")[4:8], # MDP, cMo, ncMo
    rev(RColorBrewer::brewer.pal(9, "BrBG")[c(3,2)]), # PreRPM,  RPM
    RColorBrewer::brewer.pal(9, "PRGn")[2:3], # cDC1
    RColorBrewer::brewer.pal(10, "PRGn")[8:10], # cDC2
    RColorBrewer::brewer.pal(9, "PiYG")[2], # MegP
    colorRampPalette(c(RColorBrewer::brewer.pal(9, "OrRd")[3], RColorBrewer::brewer.pal(9, "OrRd")[9]))(13)[1:4], # MEP
    colorRampPalette(c(RColorBrewer::brewer.pal(9, "OrRd")[3], RColorBrewer::brewer.pal(9, "OrRd")[9]))(13)[5:8], # ProEB
    colorRampPalette(c(RColorBrewer::brewer.pal(9, "OrRd")[3], RColorBrewer::brewer.pal(9, "OrRd")[9]))(13)[9:13] # EB
    
)

names(cell_type_fine) <- c(
    
    "MLP", 
    
    "GMP", 
    "BasoP", 
    "MastP", 
    
    "MDP",
    "cMo (1)",
    "cMo (2)",
    "ncMo (1)",
    "ncMo (2)", 

    "RPM",
    "PreRPM",
    
    "cDC1 (1)", 
    "cDC1 (2)",
    
    "cDC2 (1)",
    "cDC2 (2)",
    "cDC2 (3)",
    
    "MegP", 
    
    "MEP (1)",
    "MEP (2)", 
    "MEP (3)", 
    "MEP (4)",
    
    "ProEB (1)",  
    "ProEB (2)",
    "ProEB (3)", 
    "ProEB (4)", 
    
    "EB (1)",
    "EB (2)",
    "EB (3)", 
    "EB (4)", 
    "EB (5)"

)

#############################
### cell_type_fine_detail ###
#############################

cell_type_fine_detail <- c(
    
    RColorBrewer::brewer.pal(9, "YlOrRd")[2], # MLP
    RColorBrewer::brewer.pal(9, "BrBG")[7:9], # GMP
    RColorBrewer::brewer.pal(9, "Blues")[4:8], # MDP, cMo, ncMo
    rev(RColorBrewer::brewer.pal(9, "BrBG")[c(3,2)]), # PreRPM,  RPM
    RColorBrewer::brewer.pal(9, "PRGn")[2:3], # cDC1
    RColorBrewer::brewer.pal(10, "PRGn")[8:10], # cDC2
    RColorBrewer::brewer.pal(9, "PiYG")[2], # MegP
    colorRampPalette(c(RColorBrewer::brewer.pal(9, "OrRd")[3], RColorBrewer::brewer.pal(9, "OrRd")[9]))(13)[1:4], # MEP
    colorRampPalette(c(RColorBrewer::brewer.pal(9, "OrRd")[3], RColorBrewer::brewer.pal(9, "OrRd")[9]))(13)[5:8], # ProEB
    colorRampPalette(c(RColorBrewer::brewer.pal(9, "OrRd")[3], RColorBrewer::brewer.pal(9, "OrRd")[9]))(13)[9:13] # EB
    
)

names(cell_type_fine_detail) <- c(
    
    "MLP", 
    
    "GMP", 
    "BasoP", 
    "MastP", 
    
    "MDP",
    "cMo Ly6c hi (1)",
    "cMo Ly6c lo (2)",
    "ncMo Cd4- (1)",
    "ncMo Cd4+ (2)", 

    "RPM",
    "PreRPM",
    
    "cDC1 Cd8+ (2)", 
    "cDC1 Cd8+ prolif. (1)", 
    
    "cDC2 (1)",
    "cDC2 prolif. (2)",
    "cDC2 Ccr7+ (3)",
    
    "MegP", 
    
    "MEP (1)",
    "MEP (2)", 
    "MEP (3)", 
    "MEP (4)",
    
    "ProEB (1)",  
    "ProEB (2)",
    "ProEB (3)", 
    "ProEB (4)", 
    
    "EB (1)",
    "EB (2)",
    "EB (3)", 
    "EB (4)", 
    "EB (5)"

)
         

# Color list for R
color <- list(
    
    cellranger_class=as.list(cellranger_class), 
    tissue=as.list(tissue), 
    treatment=as.list(treatment), 
    label_main_haemosphere=as.list(label_main_haemosphere), 
    label_fine_haemosphere=as.list(label_fine_haemosphere),
    cc_phase_class=as.list(cc_phase_class), 
    cell_type_main=as.list(cell_type_main), 
    cell_type_fine=as.list(cell_type_fine), 
    cell_type_fine_detail=as.list(cell_type_fine_detail)

)