library_load <- suppressMessages(
    
    suppressWarnings(
        
        list(

            library(ggplot2), 
            library(RColorBrewer)

        )
    
    )
    
)

######################
### theme_set_text ###
######################
theme_global_set <- function(size_select=2) {
    
    size <- data.frame(
    
        size_1=c(16, 18, 20),
        size_2=c(10, 12, 16),
        size_3=c(8, 10, 12), 
        size_4=c(6, 6, 8)
        
    )

    size <- size[, size_select]
  
    theme(
        
        panel.background=element_blank(), 
        panel.spacing=unit(0.2, "lines"),

        plot.title=element_text(size=size[3], face="plain", margin=margin(t=0, r=0, b=10, l=0), color="black"), 
        plot.subtitle=element_text(size=size[2], face="plain", margin=margin(t=0, r=0, b=10, l=0), color="black"), 
    
        axis.title.y=element_text(size=size[2], face="plain", margin=margin(t=0, r=5, b=0, l=0), angle=90, color="black"), 
        axis.title.x=element_text(size=size[2], face="plain", margin=margin(t=5, r=0, b=0, l=0), color="black"),
        axis.text.y=element_text(size=size[1], face="plain", margin=margin(t=0, r=2, b=0, l=0), color="black"),
        axis.text.x=element_text(size=size[1], face="plain", margin=margin(t=2, r=0, b=0, l=0), color="black"),
      
        axis.ticks=element_line(color="black", linewidth=unit(1/2.141959, "pt")), 
        axis.line=element_line(color="black", linewidth=unit(1/2.141959, "pt")), 
      
        legend.key.height=unit(0.5, "cm"), 
        legend.key.width=unit(0.5, "cm"),
        legend.key.size=unit(0.5, "cm"),
      
        legend.title=element_text(size=size[2], face="plain", color="black"),
        legend.text=element_text(size=size[1], face="plain", color="black", hjust=0), 
        legend.key=element_rect(fill="transparent", colour="transparent"), 

        strip.text=element_text(size=size[2], margin=margin(t=2, r=2, b=2, l=2), face="plain", color="black"), 
        strip.background=element_rect(fill="white", color="black", linewidth=unit(1/2.141959, "pt"))
        
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

# FACS
facs <- c(RColorBrewer::brewer.pal(9, "BrBG")[2], RColorBrewer::brewer.pal(9, "BrBG")[8], RColorBrewer::brewer.pal(9, "BrBG")[4])
names(facs) <- c("Myeloid", "Progenitor", "Mix")

# Infection 
infection <- c("#66C2A5", "#98F5D6", "#FFA0FF")
names(infection) <- c("Baseline", "NaCl", "CpG")

# Sample group
sample_group <- c("#66c2a5", "#cd34b5", "#66c2a5", "#cd34b5", "#cd34b5", "#cd34b5", "#00634A", "#FFAC1E", "#FFAC1E", "#FFAC1E")
names(sample_group) <- c("Bl6_NaCl_D6", "Bl6_CpG_D6", "IFNAR_fl_Baseline_D0", "IFNAR_fl_CpG_D1", "IFNAR_fl_CpG_D3", "IFNAR_fl_CpG_D6", "IFNAR_fl_LysM_cre_Baseline_D0", "IFNAR_fl_LysM_cre_CpG_D1", "IFNAR_fl_LysM_cre_CpG_D3", "IFNAR_fl_LysM_cre_CpG_D6")

# Genotype 
genotype <- c("#732E8F", "#D390BF", "#FA8775")
names(genotype) <- c("Bl6", "IFNAR_fl", "IFNAR_fl_LysM_cre")

# Genotype 
dpi <- c(RColorBrewer::brewer.pal(9, "Blues")[c(4, 6, 8, 9)])
names(dpi) <- c("D0", "D1", "D3", "D6")

# Cell cycle phase 
msCC_class_RNA <- c(
    
    RColorBrewer::brewer.pal(8, "Accent")[1:3]

)

names(msCC_class_RNA) <- c(

    "G1", "S", "G2M"
    
)

# Reactome top level 
top_level_color <- c(
    
    as.character(paletteer::paletteer_d("khroma::berlin"))[seq(1, 256, length.out=14)], "#7f7f7f"

)
names(top_level_color) <- c(
    
    "Cell Cycle", "DNA Repair", "DNA Replication", "Metabolism", "Metabolism of RNA", "Metabolism of proteins", "Gene expression (Transcription)", 
    "Cellular responses to stimuli", "Signal Transduction", "Immune System", "Programmed Cell Death", "Autophagy", "Vesicle-mediated transport", "Drug ADME", 
    "Other"

)

# Cell type low
celltype_low <- c(
    
    "#a2dabd",
    "#a6a76e",
    "#fed976",
    "#ffae19",
    "#468f00",
    "#cbada4",
    "#ff9977",
    "#ed492d",
    "#7f0000",
    "#00ffff",
    "#609fe8",
    "#08306b",
    "#00a697",
    "#cc417b",
    "#ffb8ff",
    "#d34fff",
    "#9d02d7",
    "#b38cb3",
    "#00d19c",
    "#81ffd5"
    
)

names(celltype_low) <- c(

    'GMP',
    'NeuP',
    'BasoP',
    'Basophil',
    'MastP',
    'MegP',
    'MEP',
    'Proerythroblast',
    'Erythroblast',
    'cMo',
    'intMo',
    'ncMo',
    'RPM',
    'cDC1',
    'cDC2',
    'moDC', 
    'cDC mig.',
    'T cell/ILC', 
    'B cell', 
    'Plasma cell'
    
)

# Color list for R
color <- list(

    facs=facs, 
    infection=infection, 
    sample_group=sample_group, 
    genotype=genotype, 
    dpi=dpi, 
    msCC_class_RNA=msCC_class_RNA, 
    celltype_low=celltype_low

)