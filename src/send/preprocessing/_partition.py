import numpy as np
import scanpy as sc
from anndata import AnnData
import pandas as pd
from itertools import product

from collections.abc import Iterable as IterableClass
from typing import Dict, Literal, Optional, Sequence, Union, Callable, List
import scanpy as sc
import os

from ._feature_selection import binomial_deviance_selection

def partition_by_group_label(adata: AnnData,
                          group_label: str,
                          model: str = 'mouse',
                          snc_subset: Optional[str] = 'intersect',
                          highly_variable_key: Optional[str] = None,
                          layer: Optional[str] = None,
                          n_top_genes: Optional[int] = 1000,
                          batch_label: Optional[str] = None,
                          sort_genes: Optional[bool] = True,
                          dim_red: Optional[str] = None) -> Dict[str, AnnData]:
    
    data_dir = os.path.join(os.path.dirname(__file__), '..', 'data')


    data_path = os.path.join(data_dir, f'{model}_snc_genes.csv.gz')
    snc_genes = pd.read_csv(data_path, index_col=0)

    adata_snc_groups = {}

    if snc_subset == 'sasp':
        snc_genes_to_use =  snc_genes[snc_genes['SASP']].index.tolist()
    elif snc_subset == 'intersect':
        snc_genes_to_use =  snc_genes[(snc_genes['QuJi'])&(snc_genes['SenNet'])].index.tolist()
    else:
        snc_genes_to_use =  snc_genes.index.tolist() # It should be the union of all genes

    groups = adata.obs[group_label].unique().tolist() 

    for group in groups:

        adata_group = adata[adata.obs[group_label] == group].copy()

        if highly_variable_key is None:
            binomial_deviance_selection(adata=adata_group,
                                        layer=layer,
                                        deviance_key = 'binomial_deviance_' + group_label,
                                        highly_variable_key = 'highly_deviant_' + group_label,
                                        n_top_genes=n_top_genes,    
                                        batch_label=batch_label,
                                        sort_genes=sort_genes)

            snc_hvg_vars = adata_group.var_names[(adata_group.var['highly_deviant_' + group_label])\
                                            &(adata_group.var_names.isin(snc_genes_to_use))]
            
        else:

            snc_hvg_vars = adata_group.var_names[(adata_group.var[highly_variable_key])\
                                            &(adata_group.var_names.isin(snc_genes_to_use))]

        adata_group =  adata_group[:, snc_hvg_vars]
    
        if dim_red is not None:
            sc.pp.neighbors(adata_group, use_rep=dim_red)
        else: # There should be so few genes that we can just generate the kNN graph on the whole gene set
            sc.pp.neighbors(adata_group, use_rep='X')
        
        adata_snc_groups[group] = adata_group

    return adata_snc_groups
