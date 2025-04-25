import numpy as np
from scipy.sparse import csr_matrix, diags
from scipy.sparse import issparse
from numpy import log1p, log, sum as npsum
from anndata import AnnData
import pandas as pd
from itertools import product

from collections.abc import Iterable as IterableClass
from typing import Dict, Literal, Optional, Sequence, Union, Callable, List

def binomial_deviance_selection(adata: AnnData, 
                                layer: str = None,
                                deviance_key: Optional[str] = 'binomial_deviance',
                                highly_variable_key: Optional[str] = 'highly_deviant',
                                n_top_genes: int = 1000,
                                batch_key: Optional[str] = None,
                                sort_genes: Optional[bool] = True) -> None:
    """
    Python implementation of the brilliantly effective feature selection method, developed by Will Townes 
    (see Townes et al. 2019: doi.org/10.1186/s13059-019-1861-6). The idea is that we use a binomial
    deviance to quantify the variability of a gene, based on a multinomial model of UMI counts. We only
    calculate the binomial deviance, whereas the scry package developed by Townes has the option to calculate
    the Poisson deviance.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.  
    highly_vairable_key : str, optional (default: 'highly_deviant')
        The key in adata.var to store the highly variable genes.
    layer : str, optional (default: None)   
        The layer of the AnnData object to use. If None, this method won't work.
    n_top_genes : int, optional (default: 1000)
        Number of top genes to select when assigning genes as 'highly variable'
    batch_key : str, optional (default: None)
        The batch label in adata.obs. If used, we calculate the binomial deviance per batch and
        then define the binomial deviance per gene as the sum of the per-batch deviances.
    sort_genes : bool, optional (default: True)
        If True, sort genes by binomial deviance.

    Returns
    -------
    adata : AnnData
        Annotated data matrix with the binomial deviance per gene stored in adata.var['binomial_deviance'].
        The top_n_genes highly variable genes are stored in adata.var[highly_variable_key].
    """

    if layer is None:    
        counts = adata.X.copy()
    else:
        counts = adata.layers[layer].copy()

    batch_keys = adata.obs[batch_key].values if batch_key is not None else None

    # Calculate the size factors as the row sums
    size_factors = counts.sum(1).A1

    binomial_deviances = calculate_binomial_deviance_batch(counts, size_factors, batch_keys)
    adata.var[deviance_key] = np.nan_to_num(binomial_deviances, nan=0) # The NaN binomial deviances shouldn't matter

    # Set the top n genes as highly variable
    idx = adata.var[deviance_key].values.argsort()[-n_top_genes:]
    mask = np.zeros(adata.var_names.shape, dtype=bool)
    mask[idx] = True
    adata.var[highly_variable_key] = mask


def calculate_binomial_deviance_batch(counts: csr_matrix,
                                      size_factors: np.ndarray,
                                      batch_keys: Optional[np.ndarray] = None,
                                      ) -> np.ndarray:
    
    if batch_keys is not None:
        batches = np.unique(batch_keys)
        n_batches = len(np.unique(batch_keys))
        binomial_deviances_per_batch = np.zeros((counts.shape[1], n_batches))

        for i, batch in enumerate(batches):
            batch_mask = batch_keys == batch
            counts_batch = counts[batch_mask, :]
            size_factors_batch = size_factors[batch_mask]
            binomial_deviances_batch = calculate_deviance(counts_batch, size_factors_batch)
            binomial_deviances_per_batch[:, i] = binomial_deviances_batch

            binomial_deviances = binomial_deviances_per_batch.sum(1)

            return binomial_deviances
    else:
        return calculate_deviance(counts, size_factors)

def calculate_deviance(counts: csr_matrix,
                        size_factors: np.ndarray,
                        ) -> np.ndarray:
    """ Calculate the binomial deviance from the counts and size factors.
    Parameters
    ----------
    counts : csr_matrix
        The UMI counts matrix. Should be cells x genes.
    size_factors : np.ndarray
        The size factors, which, in this case, are just the total UMI counts per cell.

    Returns
    -------
    deviance: np.ndarray
        The binomial deviance values for each gene.
    """
    LP = L1P = diags(size_factors**(-1)) @ counts
    LP.data = np.log(LP.data)  # log transform nonzero elements only
    L1P.data = np.log1p(-L1P.data)  # -Inf if only a single gene nonzero in a cell
    ll_sat = (counts.multiply(LP - L1P) + diags(size_factors) @ L1P).sum(axis=0).A1
    sz_sum = size_factors.sum()
    feature_sums = counts.sum(axis=0).A1
    p = feature_sums / sz_sum
    lp = np.log(p) 
    l1p = np.log1p(-p)
    ll_null = feature_sums * (lp - l1p) + sz_sum * l1p
    deviance = 2.0 * (ll_sat - ll_null)
    return deviance