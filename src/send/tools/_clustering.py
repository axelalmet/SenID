import numpy as np
import scanpy as sc
from anndata import AnnData
import pandas as pd
from itertools import product

from collections.abc import Iterable as IterableClass
from typing import Dict, Literal, Optional, Sequence, Union, Callable, List
import scanpy as sc

from ..preprocessing import partition_by_group_label

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


# def generate_snc_clusters(adata: AnnData,
#                           snc_group_label: str = 'snc_cluster',
#                           leiden_min: float = 0.05,
#                           leiden_step: float = 0.05,
#                           leiden_max: float = 0.5):
    
