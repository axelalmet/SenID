import numpy as np
import scanpy as sc
from anndata import AnnData
import pandas as pd
from itertools import product

from collections.abc import Iterable as IterableClass
from typing import Dict, Literal, Optional, Sequence, Union, Callable, List
import scanpy as sc

from ..preprocessing import binomial_deviance_selection

def generate_snc_clusters(adata: AnnData,
                          model: str = 'mouse',
                          group_label: Optional[str] = None,
                          snc_subset: Optional[str] = 'intersect',
                          highly_variable_key: Optional[str] = None,
                          layer: Optional[str] = None,
                          n_top_genes: Optional[int] = 1000,
                          batch_label: Optional[str] = None,
                          dim_red: Optional[str] = None,
                          leiden_min: float = 0.05,
                          leiden_step: float = 0.05,
                          leiden_max: float = 0.5) -> Dict[str, AnnData]:
    
    if model not in ['human', 'mouse']:
        raise ValueError("Model must be either 'mouse' or 'human'. Working on expanding to more organisms.")
    
    if snc_subset not in ['intersect', 'union', 'sasp']:
        raise ValueError("snc_subset must be either 'intersect' or 'union', or 'sasp'.")
    
    if group_label is not None:
        print(f"Partioning data based on: {group_label}")

        adata_snc_groups = partition_by_group_label(adata,
                                                    group_label,
                                                    model,
                                                    snc_subset,
                                                    highly_variable_key,
                                                    layer,
                                                    n_top_genes,
                                                    batch_label,
                                                    dim_red)

        print(f"Generated {len(adata_snc_groups)} partitions: {adata_snc_groups}")

        print("Generating SnC-based clusters for each partition...")

def partition_by_group_label(adata: AnnData,
                          group_label: str,
                          model: str = 'mouse',
                          snc_subset: Optional[str] = 'intersect',
                          highly_variable_key: Optional[str] = 'highly_variable',
                          layer: Optional[str] = None,
                          n_top_genes: Optional[int] = 1000,
                          batch_label: Optional[str] = None,
                          sort_genes: Optional[bool] = True,
                          dim_red: Optional[str] = None) -> Dict[str, AnnData]:
    
    adata_snc_groups = {}

    snc_genes = pd.read_csv(f'data/{model}_snc_genes.csv', index_col=0)

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
            binomial_deviance_selection(adata_group,
                                        layer,
                                        n_top_genes,
                                        batch_label,
                                        sort_genes)

            snc_hvg_vars = adata_group.var_names[(adata_group.var['highly_deviant'])\
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

# def generate_snc_clusters(adata: AnnData,
#                           snc_group_label: str = 'snc_cluster',
#                           leiden_min: float = 0.05,
#                           leiden_step: float = 0.05,
#                           leiden_max: float = 0.5):
    
