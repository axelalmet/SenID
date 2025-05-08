from typing import Literal, Sequence
from anndata import AnnData
import pandas as pd
import os
from cnmf import cNMF
import numpy as np
import pingouin as pg

def subset_for_sasp_communication(adata: AnnData,
                                  send_output_key: str,
                                  pval_threshold: float = None,
                                  layer: str = 'counts',
                                  model: Literal['mouse', 'human'] = 'human') -> AnnData:


    data_dir = os.path.join(os.path.dirname(__file__), '..', 'data')

    if model not in ['human', 'mouse']:
        raise ValueError("Model must be either 'mouse' or 'human'. Working on expanding to more organisms.")
    
    data_path = os.path.join(data_dir, f'{model}_snc_genes.csv.gz')

    snc_genes = pd.read_csv(data_path, index_col=0)
    sasp_genes = snc_genes[snc_genes['SASP']].index.tolist()

    if send_output_key not in adata.uns:
        raise ValueError(f"Key '{send_output_key}' not found in adata.uns.")
    if layer not in adata.layers:
        raise ValueError(f"Layer '{layer}' not found in adata.layers.")
    
    send_output = adata.uns[send_output_key]

    if pval_threshold
        send_output = send_output[send_output['pval'] < pval_threshold]


    send_sasp_ligs = list({gene for rec in send_output['ligand'].unique() for gene in rec.split('+') if gene in sasp_genes})
    send_sasp_recs = list({gene for rec in send_output['receptor'].unique() for gene in rec.split('+') if gene in sasp_genes})

    sasp_ccc_genes = list(set(send_sasp_ligs).union(set(send_sasp_recs)))

    adata_sasp = adata[:, sasp_ccc_genes]
    adata_sasp.X = adata_sasp.layers[layer]

    adata_sasp = adata_sasp[:, adata_sasp.X.sum(0) != 0]  
    adata_sasp = adata_sasp[adata_sasp.X.sum(1) != 0]

    print(f'Removed { adata.n_obs - adata_sasp.n_obs} cells with zero SASP counts.')
    print(f'Removed {len(sasp_ccc_genes)- adata_sasp.n_vars} genes from the initial SASP gene set.')

    return adata_sasp

def run_cnmf(output_dir: str,
             adata_fname: str,
             output_name: str,
             n_modules: Sequence[int]|int,
             combine_only: bool = False,
             worker_i: int = 0,
             total_workers: int = 1,
             **kwargs) -> None:
    
    cnmf_obj = cNMF(output_dir=output_dir, name=output_name)

    if combine_only:
        cnmf_obj.combine()
        cnmf_obj.k_selection_plot()

    else:
        cnmf_obj.prepare(counts_fn=f'{output_dir}/{adata_fname}', components=n_modules, **kwargs)
        cnmf_obj.factorize(worker_i=worker_i, total_workers=total_workers)
        cnmf_obj.combine()
        cnmf_obj.k_selection_plot()

def infer_consensus_modules(adata_cnmf: AnnData,
                            output_dir: str,
                            output_name: str,
                            optimal_k: int,
                            density_threshold: float = 0.1) -> None:

    cnmf_obj = cNMF(output_dir=output_dir, name=output_name)
    cnmf_obj.consensus(k=optimal_k, density_threshold=density_threshold)
    usage, spectra_scores, spectra_tpm, top_genes = cnmf_obj.load_results(K=optimal_k, density_threshold=density_threshold)

    adata_cnmf.obsm['X_cNMF'] = usage
    adata_cnmf.varm['cNMF_spectra_scores'] = spectra_scores
    adata_cnmf.varm['cNMF_spectra_tpm'] = spectra_tpm
    adata_cnmf.uns['cNMF_top_genes'] = top_genes

# Analysis (ligand enrichment), 
def calculate_ligand_enrichment(adata: AnnData,
                                spectra_tpm_key: str = 'cNMF_spectra_tpm') -> pd.DataFrame:

    if spectra_tpm_key not in adata.varm:
        raise ValueError(f"Key '{spectra_tpm_key}' not found in adata.varm.")
    
    spectra_tpm = adata.varm[spectra_tpm_key]

#  Calculate effect sizes (log-fold change, cohen's D)
# 
