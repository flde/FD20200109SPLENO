import scanpy as sc
import pandas as pd
import numpy as np
import os


os.chdir('/research/peer/fdeckert/FD20200109SPLENO')

adata = sc.read_h5ad('data/object/qc.h5ad')
obs = pd.read_csv('data/object/components/meta.csv', index_col=0)

adata = adata[adata.obs.index.isin(obs.index.tolist())]
obs = obs.reindex(adata.obs_names)

adata.obs = obs

cell_type_fine_ery = [
    
    'MEP (1)',
    'MEP (2)', 
    'MEP (3)', 
    'MEP (4)',
    
    'ProEB (1)',  
    'ProEB (2)',
    'ProEB (3)', 
    'ProEB (4)', 
    
    'EB (1)',
    'EB (2)',
    'EB (3)', 
    'EB (4)', 
    'EB (5)'
    
    
]

adata_1 = adata[adata.obs['cell_type_fine'].isin(cell_type_fine_ery)].copy()
adata_1.obs['cell_type_fine'] = pd.Series(adata_1.obs['cell_type_fine'], dtype='category')

sc.pp.normalize_total(adata_1)
sc.pp.log1p(adata_1)
sc.pp.scale(adata_1)
sc.tl.pca(adata_1, svd_solver='arpack')
sc.pp.neighbors(adata_1, n_neighbors=20, use_rep='X_pca')

sc.tl.diffmap(adata_1, n_comps=15)

adata_1.uns['iroot'] = np.flatnonzero(adata_1.obs['cell_type_fine']=='MEP (1)')[0]
sc.tl.dpt(adata_1, n_branchings=0, n_dcs=10)

pd.DataFrame(adata_1.obsm['X_diffmap'], index=adata_1.obs.index).to_csv('data/object/components/diffusion_map.csv')
pd.DataFrame(adata_1.obs['dpt_pseudotime'], index=adata_1.obs.index).to_csv('data/object/components/diffusion_pseudotime.csv')