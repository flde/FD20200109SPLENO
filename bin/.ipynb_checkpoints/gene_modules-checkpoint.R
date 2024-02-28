########################################################
### Progenitor marker genes from Dahlin et al., 2018 ###
########################################################
genes_hsc <- c("Procr") # Not expressed: Fgd5 and Hoxb5
genes_ly <- c("Dntt", "Flt3") # Not expressed: Rag2 
genes_meg <- c("Pf4", "Itga2b", "Gp1bb")
genes_ery <- c("Klf1", "Epor", "Gata1")
genes_mop <- c("Irf8", "Csf1r", "Ly86") 
genes_neu <- c("Elane", "Gfi1")
genes_eo <- c("Prg2", "Prg3") 
genes_baso <- c("Prss34", "Mcpt8")
genes_mast <- c("Gzmb", "Cma1")

# Additional marker 
genes_hsc <- c(genes_hsc, "Cd34") # Wilson et al. 2015 (Ly6a)
genes_ly <- c(genes_ly) 
genes_meg <- c(genes_meg)
genes_ery <- c(genes_ery, "Hba-a2") 
genes_mop <- c(genes_mop, "Spi1", "Irf8", "Klf4") # PU.1 (Spi1)
genes_mo <- c("Ms4a3", "Cx3cr1", "Flt3", "Csf1r", "Ccr2", "Adgre1") # Ms4a3 https://www.sciencedirect.com/science/article/pii/S0092867419308979?via%3Dihub
genes_neu <- c(genes_neu)
genes_eo <- c(genes_eo)
genes_baso <- c(genes_baso)
genes_mast <- c(genes_mast)
genes_rpm <- c("Spic", "Vcam1")
genes_dc1 <- c("Cd8a", "Xcr1", "Itgae") # https://www.sciencedirect.com/science/article/pii/S009286741931116X?via%3Dihub
genes_dc2 <- c("Zbtb46", "Cd4", "Dtx1", "Esam", "Tbx21") # Tbet:Tbx21. The Tbet- cDC2 subset should express Clec12a but all myeloid express that. Zbtb46 Should disting form Mo/Mac


# Combine results 
genes_marker <- data.frame(
    
    cell_type=factor(
        
        c(
            rep("HSC", length(genes_hsc)), 
            rep("Meg", length(genes_meg)), 
            rep("Ery", length(genes_ery)), 
            rep("Neu", length(genes_neu)), 
            rep("Eo", length(genes_eo)), 
            rep("Mast", length(genes_mast)), 
            rep("Baso", length(genes_baso)), 
            rep("GMP/MDP", length(genes_mo)), 
            rep("RPM", length(genes_rpm)), 
            rep("cDC1", length(genes_dc1)), 
            rep("cDC2", length(genes_dc2))
    ), 
        
        levels=c("HSC", "Meg", "Ery", "Neu", "Eo", "Mast", "Baso", "GMP/MDP", "RPM", "cDC1", "cDC2")
    
    ), 
    
    genes=c(genes_hsc, genes_meg, genes_ery, genes_neu, genes_eo, genes_mast, genes_baso, genes_mo, genes_rpm, genes_dc1, genes_dc2)


)
###################################
### Figure markers for dot plot ###
###################################
genes_hsc <- c("Procr", "Cd34") # Not expressed: Fgd5 and Hoxb5
genes_meg <- c("Pf4", "Itga2b")
genes_ery <- c("Klf1", "Epor", "Hba-a2")
genes_neu <- c("Elane", "Gfi1")
genes_eo <- c("Prg2", "Prg3") 
genes_mast <- c("Gzmb", "Cma1")
genes_baso <- c("Prss34", "Mcpt8")
genes_gmpdmp <- c("Ms4a3", "Cx3cr1", "Flt3", "Ly6c2") # Ms4a3 https://www.sciencedirect.com/science/article/pii/S0092867419308979?via%3Dihub
genes_rpm <- c("Spic", "Adgre1", "Cd163")
genes_dc1 <- c("Cd8a", "Xcr1", "Itgae") # https://www.sciencedirect.com/science/article/pii/S009286741931116X?via%3Dihub
genes_dc2 <- c("Cd4", "Esam", "Zbtb46") # Tbet:Tbx21. The Tbet- cDC2 subset should express Clec12a but all myeloid express that. Zbtb46 Should disting form Mo/Mac


# Combine results 
genes_marker_figure <- data.frame(
    
    cell_type=factor(
        
        c(
            rep("HSC", length(genes_hsc)), 
            rep("MegP", length(genes_meg)), 
            rep("Ery", length(genes_ery)), 
            rep("NeuP", length(genes_neu)), 
            rep("EoP", length(genes_eo)), 
            rep("MastP", length(genes_mast)), 
            rep("BasoP", length(genes_baso)), 
            rep("GMP/MDP", length(genes_gmpdmp)), 
            rep("RPM", length(genes_rpm)), 
            rep("cDC1", length(genes_dc1)), 
            rep("cDC2", length(genes_dc2))
    ), 
        
        levels=c("HSC", "MegP", "Ery", "NeuP", "EoP", "MastP", "BasoP", "GMP/MDP", "RPM", "cDC1", "cDC2")
    
    ), 
    
    genes=c(genes_hsc, genes_meg, genes_ery, genes_neu, genes_eo, genes_mast, genes_baso, genes_gmpdmp, genes_rpm, genes_dc1, genes_dc2)


)

#####################################
### HSCs from Wilson et al., 2015 ###
#####################################

genes_hsc_wilson <- c(
    
    "Slamf1", 
    "Procr", 
    "Pdzk1ip1", 
    "Ltb", 
    "Mllt3", 
    "Ifitm1", 
    "Gimap1", 
    "Gimap6", 
    "Limd2", 
    "Trim47", 
    "Neil2", 
    "Vwf", 
    "Pde1b", 
    "Neo1", 
    "Sqrdl", 
    "Sult1a1", 
    "Cd82", 
    "Ramp2", 
    "Ubl3", 
    "Ly6a", 
    "Cdkn1c", 
    "Fgfr3", 
    "Cldn10", 
    "Ptpn14", 
    "Mettl7a1", 
    "Smtnl1", 
    "Ctsf", 
    "Gstm1", 
    "Sox18", 
    "Fads3"
    
)


##################################
### REACTOME pathway selection ###
##################################

reactome_filter_eb <- c(

    "REACTOME_SCAVENGING_OF_HEME_FROM_PLASMA", 
#     "REACTOME_ERYTHROCYTES_TAKE_UP_OXYGEN_AND_RELEASE_CARBON_DIOXIDE", 
#     "REACTOME_ERYTHROCYTES_TAKE_UP_CARBON_DIOXIDE_AND_RELEASE_OXYGEN", 
    "REACTOME_HEME_SIGNALING", 
    "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING", 
    "REACTOME_TP53_REGULATES_METABOLIC_GENES", 
    "REACTOME_HEME_BIOSYNTHESIS", 
    "REACTOME_DDX58_IFIH1_MEDIATED_INDUCTION_OF_INTERFERON_ALPHA_BETA", 
    "REACTOME_INTERFERON_GAMMA_SIGNALING", 
    "REACTOME_MITOTIC_G2_G2_M_PHASES", 
    "REACTOME_HEDGEHOG_ON_STATE", 
    "REACTOME_HEDGEHOG_OFF_STATE", 
    "REACTOME_TCF_DEPENDENT_SIGNALING_IN_RESPONSE_TO_WNT", 
    "REACTOME_SIGNALING_BY_HEDGEHOG", 
    "REACTOME_IRON_UPTAKE_AND_TRANSPORT", 
    "REACTOME_G2_M_CHECKPOINTS"
    
)