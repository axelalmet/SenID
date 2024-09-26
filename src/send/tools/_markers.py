import numpy as np
import scanpy as sc
from anndata import AnnData
import pandas as pd
from itertools import product
from scipy.sparse import vstack
from scipy.stats import mannwhitneyu

from collections.abc import Iterable as IterableClass
from typing import Dict, Literal, Optional, Sequence, Union, Callable, List

def score_marker_genes(adata: AnnData,
                        group_label: str,
                        comparison: str = 'pairwise') -> pd.DataFrame:
    """
    """

    if comparison not in ['pairwise', 'one_vs_rest']:
        raise ValueError('Comparison must be either "pairwise" or "one_vs_rest".')
    
    gene_vars = adata.var_names

    expr_norm = np.expm1(adata.X)
    expr_norm_squared = expr_norm.copy()
    expr_norm_squared.data **= 2

    groups = adata.obs[group_label].unique().tolist()

    lfc_results = {group: {} for group in groups}

    if comparison == 'pairwise':

        group_pairs = [pair for pair in product(groups, repeat=2) if pair[0] != pair[1]]

        for pair in group_pairs:
            group1, group2 = pair
            group1_indices = np.ravel(adata.obs[group_label] == group1)
            group2_indices = np.ravel(adata.obs[group_label] == group2)
            

            results_df = _calculate_effect_sizes(expr_norm,
                            expr_norm_squared,
                            group1_indices,
                            group2_indices,
                            gene_vars)
            
            results_df['comparison'] =  f'{group1}_vs_{group2}'

            lfc_results[group1][group2] = results_df

        # Aggregate the results per group


    else:

        for group in groups:
            group1_indices = np.ravel(adata.obs[group_label] == group)
            group2_indices = ~group1_indices

            results_df = _calculate_effect_sizes(expr_norm,
                            expr_norm_squared,
                            group1_indices,
                            group2_indices,
                            gene_vars)
            
            results_df['comparison'] =  f'{group1}_vs_rest'
            
            lfc_results[group1] = results_df

def _calculate_effect_sizes(expr_norm: np.ndarray,
                            expr_norm_squared: np.ndarray,
                            group1_indices: np.ndarray,
                            group2_indices: np.ndarray,
                            gene_vars: Union[List[str], pd.Index]) -> pd.DataFrame:
    """
    """
    n1 = group1_indices.sum()
    n2 = group2_indices.sum()

    group1_expr = expr_norm[group1_indices]
    group2_expr = expr_norm[group2_indices]

    group1_detected = (group1_expr > 0).sum(0) / n1
    group2_detected = (group2_expr > 0).sum(0) / n2

    pooled_expr = vstack([group1_expr, group2_expr])

    group1_mean = group1_expr.mean(axis=0).A1
    group2_mean = group2_expr.mean(axis=0).A1

    group1_var = (expr_norm_squared[group1_indices].mean(0) - np.square(group1_expr.mean(0))).A1
    group2_var = (expr_norm_squared[group2_indices].mean(0) - np.square(group2_expr.mean(0))).A1

    lfc = np.log2(group1_mean + 1e-10) - np.log2(group2_mean + 1e-10)
    lfc_cohen = lfc / np.sqrt( 0.5 * (group1_var + group2_var) )

    lfc_detected = np.log2(group1_detected.A1 + 1e-10) - np.log2(group2_detected.A1 + 1e-10)

    # Calculate the mann-whitney U statistic so we can calculate the AUC
    u_statistic = mannwhitneyu(group1_expr.toarray(), group2_expr.toarray(), axis=0).statistic
    auc = u_statistic/(group1_expr.shape[0] * group2_expr.shape[0])

    # Let us rank 
    results_df = pd.DataFrame({'logfoldchanges': lfc, 'logfoldchanges_cohen': lfc_cohen,
                                'logfoldchanges_detected': lfc_detected, 'auc': auc})
    
    results_df.index = pd.Index(gene_vars)
    
    return results_df
