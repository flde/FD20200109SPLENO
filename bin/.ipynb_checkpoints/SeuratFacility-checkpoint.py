#################
### dir2adata ###
#################

def dir2adata(path, assay, slot):
    
    import scanpy as sc
    import pandas as pd
    
    import os 
    
    os.chdir('/research/peer/fdeckert/FD20200109SPLENO')
    
    # Create anndata from count mtx wit var_names 
    adata = sc.read_mtx(path+'assay/'+assay+'/'+slot+'/'+'matrix.mtx')
    adata.var_names = pd.read_table(path+'assay/'+assay+'/'+slot+'/'+'genes.csv', index_col=0).index
    
    # Add obs from meta_data
    adata.obs = pd.read_csv(path+'meta/meta.csv', index_col=0)
    
    # Set layers
    adata.layers['counts'] = adata.X
    
    return(adata)

################
### BuildDir ###
################

def BuildDir(path, overwrite=False): 
    
    import os
    import sys
    import shutil
    
    if os.path.exists(path) and overwrite==False:
        
        print('Directory already exists and overwrite is set to False')
    
    else: 
        
        print('Creating output directory'+path)
        
        if os.path.exists(path): shutil.rmtree(path)
        os.makedirs(path)
        
#####################
### WriteMetaData ###
#####################
        
def WriteMetaData(adata, path, overwrite=False, file='meta.csv'): 
    
    import os
    import pandas as pd
    
    if os.path.exists(path+'meta/'+file) and overwrite==False:
        
        print('Writing meta data: The meta data file already exists and overwrite is set to False')
        
    else: 
        
        # Create meta data directory
        if os.path.exists(path+'meta/')==False: 
            os.makedirs(path+'meta/')
        
        # Write meta data 
        print('Writing meta data: '+path+'meta/'+file)
        adata.obs.to_csv(path+'meta/'+file)
        
######################
### WriteAssayData ###
######################

def WriteAssayData(adata, path, assay, overwrite=False): 
    
    import os
    import scipy as sc
    
    print('Writing assay data for the slots:')

    for slot in adata.layers:
        
        if os.path.exists(path+'assay/'+assay+'/'+slot) and overwrite==False: 
            
            print('... '+assay+' '+slot+' already exists and overwrite is set to False')
            
        else: 
            
            print('... ', assay, slot)
            if os.path.exists(path+'assay/'+assay+'/'+slot)==False: 
                os.makedirs(path+'assay/'+assay+'/'+slot)
            
            sc.io.mmwrite(path+'assay/'+assay+'/'+slot+'/matrix.mtx', sc.sparse.csc_matrix(adata.layers[slot]))
            adata.var_names.to_series().to_csv(path+'assay/'+assay+'/'+slot+'/genes.csv')
            adata.obs_names.to_series().to_csv(path+'assay/'+assay+'/'+slot+'/cellid.csv')
            
##########################
### WriteReductionData ###
##########################

def WriteReductionData(adata, path, overwrite=False): 
    
    import os
    import numpy as np
    import pandas as pd
    import scipy as sc
    
    print('Writing data for the reduction:')

    for reductions in adata.obsm:
        
        if os.path.exists(path+'reductions/'+reductions) and overwrite==False: 
            
            print('... '+reductions+' already exists and overwrite is set to False')
            
        else: 
            
            print('... ', reductions)
            if os.path.exists(path+'reductions/'+reductions)==False: 
                os.makedirs(path+'reductions/'+reductions)

            pd.DataFrame(adata.obsm[reductions], index=adata.obs_names.to_series()).to_csv(path+'reductions/'+reductions+'/reduction.csv')

#################
### adata2dir ###
#################

def adata2dir(adata, path, assay, slot, build_dir=False, overwrite=False): 
    
    import os
    
    os.chdir('/research/peer/fdeckert/FD20200109SPLENO')
    
    if(build_dir==True):
        
        BuildDir(path, overwrite)

    # Write meta data 
    WriteMetaData(adata, path, overwrite)

    # Write assay data 
    WriteAssayData(adata, path, assay, overwrite)
    
    # Write reduction data
    WriteReductionData(adata, path, overwrite)
    












