##################
### cell types ###
##################
cell_type <- data.frame(

    ident=c(
        
        "seurat_clusters_eb_13", 
        "seurat_clusters_eb_1", 
        "seurat_clusters_eb_12", 
        "seurat_clusters_eb_11", 
        "seurat_clusters_eb_9", 
        "seurat_clusters_eb_0", 
        "seurat_clusters_eb_5", 
        "seurat_clusters_eb_3", 
        "seurat_clusters_eb_7", 
        "varid_clusters_mpp_8", 
        "varid_clusters_mpp_7", 
        "varid_clusters_mpp_6", 
        "varid_clusters_mpp_4", 
        "varid_clusters_mpp_1", 
        "varid_clusters_mpp_3", 
        "varid_clusters_mpp_2", 
        "varid_clusters_mpp_5", 
        "varid_clusters_mpp_9", 
        "varid_clusters_mpp_10", 
        "varid_clusters_m_1", 
        "varid_clusters_m_7", 
        "varid_clusters_m_6", 
        "varid_clusters_m_3", 
        "varid_clusters_m_8", 
        "varid_clusters_m_5", 
        "varid_clusters_m_4", 
        "varid_clusters_m_2", 
        "varid_clusters_m_10", 
        "varid_clusters_m_11", 
        "varid_clusters_m_9"
    
    ), 
    
    cell_type_main=c(
        
        "ProEB",
        "EB",
        "ProEB",
        "ProEB",
        "EB",
        "EB",
        "ProEB",
        "EB",
        "EB",
        "MDP",
        "MLP",
        "GMP",
        "GMP",
        "MEP",
        "MEP",
        "MEP",
        "MegP",
        "MEP",
        "GMP",
        "ncMo",
        "ncMo",
        "RPM",
        "cMo",
        "cDC2",
        "cDC1",
        "cDC1",
        "cDC2",
        "RPM",
        "cMo",
        "cDC2" 
    
    ), 
    
    cell_type_fine=c(
        
        "ProEB (1)",
        "EB (5)",
        "ProEB (3)",
        "ProEB (4)",
        "EB (1)",
        "EB (4)",
        "ProEB (2)",
        "EB (3)",
        "EB (2)",
        "MDP",
        "MLP",
        "BasoP",
        "GMP",
        "MEP (2)",
        "MEP (3)",
        "MEP (4)",
        "MegP",
        "MEP (1)",
        "MastP",
        "ncMo (1)",
        "ncMo (2)",
        "RPM",
        "cMo (2)",
        "cDC2 (3)",
        "cDC1 (2)",
        "cDC1 (1)",
        "cDC2 (1)",
        "PreRPM",
        "cMo (1)",
        "cDC2 (2)"

    ), 
    
    cell_type_fine_detail=c(
        
        "ProEB (1)",
        "EB (5)",
        "ProEB (3)",
        "ProEB (4)",
        "EB (1)",
        "EB (4)",
        "ProEB (2)",
        "EB (3)",
        "EB (2)",
        "MDP",
        "MLP",
        "BasoP",
        "GMP",
        "MEP (2)",
        "MEP (3)",
        "MEP (4)",
        "MegP",
        "MEP (1)",
        "MastP",
        "ncMo Cd4- (1)",
        "ncMo Cd4+ (2)",
        "RPM",
        "cMo Ly6c lo (2)",
        "cDC2 Ccr7+ (3)",
        "cDC1 Cd8+ (2)",
        "cDC1 Cd8+ prolif. (1)",
        "cDC2 (1)",
        "PreRPM",
        "cMo Ly6c hi (1)",
        "cDC2 prolif. (2)"

    )

)

cell_type_fine_order <- c(
    
    "cDC1 (1)",
    "cDC1 (2)",
    "cDC2 (1)",
    "cDC2 (2)",
    "cDC2 (3)",
    "ncMo (1)",
    "ncMo (2)",
    "RPM",
    "PreRPM",
    "cMo (1)",
    "cMo (2)",
    "MDP",
    "BasoP",
    "GMP",
    "MastP",
    "MLP",
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

cell_type_fine_detail_order <- c(
    
    "cDC1 Cd8+ (2)",
    "cDC1 Cd8+ prolif. (1)",
    "cDC2 (1)",
    "cDC2 prolif. (2)", 
    "cDC2 Ccr7+ (3)",
    "PreRPM",
    "RPM",
    "ncMo Cd4- (1)",
    "ncMo Cd4+ (2)", 
    "cMo Ly6c hi (1)",
    "cMo Ly6c lo (2)",
    "MDP",
    "BasoP",
    "GMP",
    "MLP",
    "MastP",
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

cell_type_main_order <- c(
    
    "cDC1",
    "cDC1",
    "cDC2",
    "cDC2",
    "cDC2",
    "PreRPM",
    "RPM",
    "ncMo",
    "ncMo",
    "cMo",
    "cMo",
    "MDP",
    "GMP",
    "GMP",
    "MLP",
    "GMP",
    "MegP",
    "MEP",
    "MEP",
    "MEP",
    "MEP",
    "ProEB",
    "ProEB",
    "ProEB",
    "ProEB",
    "EB",
    "EB",
    "EB",
    "EB",
    "EB"

) 

cell_type <- cell_type[match(cell_type_fine_order, cell_type$cell_type_fine),]

###############
### Subsets ###
###############
cell_type_prog <- c(

    "MDP",
    "BasoP",
    "GMP",
    "MLP",
    "MastP",
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