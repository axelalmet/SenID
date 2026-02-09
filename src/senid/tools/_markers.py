from __future__ import annotations

import numpy as np
import scanpy as sc
from anndata import AnnData
import pandas as pd
from pandas import DataFrame
from itertools import product
from scipy.sparse import vstack, csr_matrix
from scipy.stats import mannwhitneyu

from collections.abc import Iterable as IterableClass
from typing import Dict, Literal, Optional, Sequence, Union, Callable, List

def score_marker_genes(adata: AnnData,
                        group_label: str,
                        output_key: str = 'SenID_marker_scores',
                        comparison: str = 'pairwise',
                        return_df: bool = False) -> DataFrame | None:
    """
    """

    if comparison not in ['pairwise', 'one_vs_rest']:
        raise ValueError('Comparison must be either "pairwise" or "one_vs_rest".')
    
    gene_vars = adata.var_names

    expr = adata.X.copy()
    expr_norm = np.expm1(adata.X)
    expr_squared = expr.copy()
    expr_squared.data **= 2

    groups = adata.obs[group_label].unique().tolist()

    lfc_results = {group: {} for group in groups}
    lfc_results_aggregated = {group: {} for group in groups}

    effect_sizes = ['logfoldchanges', 'logfoldchanges_cohen', 'logfoldchanges_detected', 'auc']
    if comparison == 'pairwise':

        group_pairs = [pair for pair in product(groups, repeat=2) if pair[0] != pair[1]]

        for pair in group_pairs:
            group1, group2 = pair
            group1_indices = np.ravel(adata.obs[group_label] == group1)
            group2_indices = np.ravel(adata.obs[group_label] == group2)
            
            results_df = _calculate_effect_sizes(expr,
                                                 expr_norm,
                                                expr_squared,
                                                group1_indices,
                                                group2_indices,
                                                gene_vars)
            results_df.fillna(value=0, inplace=True)
            
            results_df['comparison'] =  f'{group1}_vs_{group2}'

            lfc_results[group1][group2] = results_df

        # Aggregate the results per group
        for group1 in lfc_results:
            # Concatenate all DataFrames in the sub-dictionary along the columns
            combined_effect_sizes = {eff_size: pd.concat([lfc_results[group1][group2][eff_size] for group2 \
                                                          in lfc_results[group1]], axis=1)\
                                                for eff_size in effect_sizes}
            
            # Use agg to calculate mean, median, min, and max
            summary_df = {eff_size: combined_effect_sizes[eff_size].agg(['mean', 'median', 'min', 'max'], axis=1)\
                                    for eff_size in effect_sizes}
            
            # Store the result in the summary_results dictionary
            group1_df = pd.DataFrame()
            group1_df['gene'] = gene_vars

            for eff_size in effect_sizes:
                group1_df['mean_' + eff_size] = summary_df[eff_size]['mean']
                group1_df['min_' + eff_size] = summary_df[eff_size]['min']
                group1_df['median_' + eff_size] = summary_df[eff_size]['median']
                group1_df['max_' + eff_size] = summary_df[eff_size]['max']

            group1_df['group'] = group1
            group1_df['comparison'] = 'pairwise'
            lfc_results_aggregated[group1] = group1_df

        marker_df = pd.concat([lfc_results_aggregated[group] for group in groups], axis=0)

    # Else we're looking at "one vs rest"
    else:

        for group in groups:
            group1_indices = np.ravel(adata.obs[group_label] == group)
            group2_indices = ~group1_indices

            results_df = _calculate_effect_sizes(expr,
                                                 expr_norm,
                                                expr_squared,
                                                group1_indices,
                                                group2_indices,
                                                gene_vars)
            results_df.fillna(value=0, inplace=True)

            results_df['group'] =  group1
            results_df['comparison'] =  'one_vs_rest'
            
            lfc_results[group1] = results_df

            marker_df = pd.concat([lfc_results[group] for group in groups], axis=0)
            marker_df.rename(columns={'logfoldchanges': 'mean_logfoldchanges',
                                      'logfoldchanges_cohen':'mean_logfoldchanges_cohen',
                                      'logfoldchanges_detected': 'mean_logfoldchanges_detected',
                                      'auc': 'mean_auc'}, inplace=True)
            
            for eff_size in effect_sizes:
                marker_df['min_' + eff_size] = marker_df['mean_' + eff_size].values
                marker_df['median_' + eff_size] = marker_df['median_' + eff_size].values
                marker_df['max_' + eff_size] = marker_df['max_' + eff_size].values


    if return_df:
        return marker_df
    else:
        adata.uns[output_key] = marker_df

def _calculate_effect_sizes(expr: np.ndarray,
                            expr_norm: np.ndarray,
                            expr_squared: np.ndarray,
                            group1_indices: np.ndarray,
                            group2_indices: np.ndarray,
                            gene_vars: Union[List[str], pd.Index]) -> pd.DataFrame:
    """
    """
    n1 = group1_indices.sum()
    n2 = group2_indices.sum()

    group1_expr = expr[group1_indices]
    group2_expr = expr[group2_indices]

    group1_expr_norm = expr_norm[group1_indices]
    group2_expr_norm = expr_norm[group2_indices]

    group1_detected = (group1_expr > 0).sum(0) / n1
    group2_detected = (group2_expr > 0).sum(0) / n2

    group1_mean = group1_expr.mean(axis=0).A1
    group2_mean = group2_expr.mean(axis=0).A1

    group1_norm_mean = group1_expr_norm.mean(axis=0).A1
    group2_norm_mean = group2_expr_norm.mean(axis=0).A1

    group1_var = (expr_squared[group1_indices].mean(0) - np.square(group1_expr.mean(0))).A1
    group2_var = (expr_squared[group2_indices].mean(0) - np.square(group2_expr.mean(0))).A1

    lfc = np.log2(group1_norm_mean + 1e-12) - np.log2(group2_norm_mean + 1e-12)
    lfc = np.nan_to_num(lfc, nan=0)# If there's a nan it means we had a 'log(0) - log(0) case

    log_means_diff = group1_mean - group2_mean
    lfc_cohen = log_means_diff / ( 0.5 * (group1_var + group2_var) )** 0.5

    lfc_detected = np.log2(group1_detected.A1 + 1e-12) - np.log2(group2_detected.A1 + 1e-12)

    # Calculate the mann-whitney U statistic so we can calculate the AUC
    u_statistic = mannwhitneyu(group1_expr_norm.toarray(), group2_expr_norm.toarray(), axis=0).statistic
    auc = u_statistic/(group1_expr_norm.shape[0] * group2_expr_norm.shape[0])
    
    results_df = pd.DataFrame({'gene': gene_vars, 'logfoldchanges': lfc, 'logfoldchanges_cohen': lfc_cohen,
                                'logfoldchanges_detected': lfc_detected, 'auc': auc})
        
    return results_df
