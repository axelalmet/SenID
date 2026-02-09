
from typing import Literal, Dict
# import tensorflow as tf  # For building and training the machine learning model
import pandas as pd  # For data manipulation and analysis
from anndata import AnnData
from scanpy import pp  # For preprocessing single-cell data
import scanpy as sc
import os
from os import path  # For pathname manipulations
# from mNSF import MoranI # For calculating spatial dependeicy of each factor within each sample, using Moran's I
from itertools import product
import numpy as np
import warnings

# Ignore all DeprecationWarnings
warnings.filterwarnings("ignore", category=DeprecationWarning)


def subset_for_spatial_sasp_communication(adata: AnnData | Dict[str, AnnData],
                                   layer: str = 'counts',
                                   db_name: str = None,
                                   model: Literal['mouse', 'human'] = 'human') -> AnnData:
    """
    Subset the AnnData object for spatial SASP communication analysis.
    """

    data_dir = os.path.join(os.path.dirname(__file__), '..', 'data')

    if model not in ['human', 'mouse']:
        raise ValueError("Model must be either 'mouse' or 'human'. Working on expanding to more organisms.")
    
    data_path = os.path.join(data_dir, f'{model}_snc_genes.csv.gz')

    snc_genes = pd.read_csv(data_path, index_col=0)
    sasp_genes = snc_genes[snc_genes['SASP']].index.tolist()
    
    if isinstance(adata, dict):
        inferred_sasp_ccc_genes = []

        for sample in adata:
            adata_sample = adata[sample]
            inferred_interactions = [col.strip('s-') for col in adata_sample.obsm[f'commot-{db_name}-sum-sender'].columns if col.count('-') > 1 and 'total' not in col]

    else:
        inferred_interactions = [col.strip('s-') for col in adata.obsm[f'commot-{db_name}-sum-sender'].columns if col.count('-') > 1 and 'total' not in col]
