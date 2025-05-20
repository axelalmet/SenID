from anndata import AnnData
from scipy.sparse import csr_matrix
import numpy as np
from itertools import product
import pandas as pd

def _determine_bursty_genes(U: csr_matrix, S: csr_matrix, var_t=1.5, u_min=0.02, s_min=0.02):
    '''
    Determine which genes show bursty (overdispersed) dynamics based on
    thresholds for mean expression and dispersion.
    '''

    def sparse_variance(X):
        # Efficient sparse variance: E[x^2] - (E[x])^2
        X_sq = X.copy()
        X_sq.data **= 2
        mean = X.mean(axis=0).A1  # .A1 flattens sparse matrix to 1D
        mean_sq = X_sq.mean(axis=0).A1
        return mean_sq - mean**2

    # Compute mean and variance
    U_mean = U.mean(axis=0).A1
    S_mean = S.mean(axis=0).A1
    U_var = sparse_variance(U)
    S_var = sparse_variance(S)

    # Dispersion
    U_disp = (U_var - U_mean) / (U_mean ** 2 + 1e-12)
    S_disp = (S_var - S_mean) / (S_mean ** 2 + 1e-12)

    # Log ratio filter (with added epsilon to avoid division-by-zero)
    log_ratio = np.abs(np.log((S_mean + 1e-12) / (U_mean + 1e-12)))

    # Combined boolean mask
    fitted_mask = (
        (U_mean > u_min) &
        (S_mean > s_min) &
        (U_disp > var_t) &
        (S_disp > var_t) &
        (log_ratio < 4)
    )

    return fitted_mask

def _filter_genes(adata: AnnData,
                        spliced_layer: str = 'spliced',
                        unspliced_layer: str = 'unspliced',
                        **kwargs) -> AnnData:
    
    if spliced_layer not in adata.layers:
        raise ValueError(f"Layer '{spliced_layer}' not found in adata.layers.")
    if unspliced_layer not in adata.layers:
        raise ValueError(f"Layer '{unspliced_layer}' not found in adata.layers.")
    
    if not isinstance(adata.layers[spliced_layer], csr_matrix):
        raise TypeError(f"Layer '{spliced_layer}' must be a sparse matrix.")
    if not isinstance(adata.layers[unspliced_layer], csr_matrix):
        raise TypeError(f"Layer '{unspliced_layer}' must be a sparse matrix.")
    
    S = adata.layers[spliced_layer]
    U = adata.layers[unspliced_layer]

    fitted_idx = _determine_bursty_genes(U, S, **kwargs)
    print('No. all genes that pass thresh: ', np.sum(fitted_idx))

    g_names_toUse = adata.var_names[fitted_idx]

    adata_monod = AnnData(X=S)
    adata_monod.layers['unspliced'] = U
    adata_monod.layers['spliced'] = S
    adata_monod.obs = adata.obs
    adata_monod.var = adata.var
    adata_monod.obs.index.name = 'barcode'
    adata_monod.var.index.name = 'gene_name'

    adata_monod = adata_monod[:, g_names_toUse]

    return adata_monod

def _determine_comparison_groups(adata: AnnData,
                                 partition_key: str,
                                 comparison_key: str,
                                 min_cells: int = 100) -> list:
    """
    Determine the groups to consider for comparison based on the number of cells.
    
    Parameters:
    - adata: AnnData object containing the data.
    - partition_key: Column name in adata.obs that contains partition information.
    - comparison_key: Column name in adata.obs that contains comparison information.
    - min_cells: Minimum number of cells required to consider a group.
    
    Returns:
    - List of groups to consider for comparison.
    """
    
    cell_groups_per_comparison = pd.crosstab(adata.obs[partition_key], adata.obs[comparison_key])

    cell_groups_to_consider = cell_groups_per_comparison[cell_groups_per_comparison.min(1) > min_cells].index.tolist()

    return cell_groups_to_consider

def generate_loom_objects(adata: AnnData,
                          dataset_key: str,
                          partition_key: str, 
                          comparison_key: str,
                          str_replace: str = '-',
                          output_directory: str = '.',
                          min_cells: int = 100,
                          return_groups: bool = False,
                          **kwargs):
    """
    Generate loom files for each phenotype in the AnnData object.
    
    Parameters:
    - adata: AnnData object containing the data.
    - pheno_col: Column name in adata.obs that contains phenotype information.
    - data_directory: Directory to save the loom files.
    """

    strs_to_replace = ['/', ' ', ',', '.']
    str_replacements = str.maketrans({str: str_replace for str in strs_to_replace})
    adata_monod = _filter_genes(adata, **kwargs)

    cell_groups = adata.obs[partition_key].unique()
    comparisons = adata.obs[comparison_key].unique()

    cell_groups_to_consider = _determine_comparison_groups(adata_monod,
                                                           partition_key,
                                                           comparison_key,
                                                           min_cells=min_cells)

    print(f'Cell groups that passed min_cells threshold: {", ".join(sorted(cell_groups_to_consider))}')
    for group, part in product(cell_groups_to_consider, comparisons):
        
        group_label = group.translate(str_replacements)
        part_label = part.translate(str_replacements)

        adata_monod_group = adata_monod[(adata_monod.obs[partition_key] == group)&(adata_monod.obs[comparison_key] == part), :].copy()

        adata_monod_group.write_loom(f'{output_directory}/{dataset_key}_{group_label}_{part_label}.loom')
        print(f'Saved loom file for {group} at {output_directory}/{dataset_key}_{group_label}_{part_label}.loom')

    if return_groups:
        return cell_groups_to_consider