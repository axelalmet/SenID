from typing import Literal, Sequence, Optional, List
from anndata import AnnData
import pandas as pd
import os
from cnmf import cNMF
import numpy as np
import pingouin as pg

def constrain_senchat_interactions(adata:  AnnData,
                                   senchat_output_key: str,
                                   implausible_interactions: dict = None,
                                   pathway_types: Optional[List[str]] = None) -> None:

    if senchat_output_key not in adata.uns:
        raise ValueError(f"Key '{senchat_output_key}' not found in adata.uns. Please run SenChat first.")

    if implausible_interactions is None and pathway_types is None:
        raise ValueError("At least one of 'implausible_interactions' or 'pathway_types' must be provided.")
    
    allowed_pathway_types = ['Secreted Signaling', 'Cell-Cell Contact', 'ECM-Receptor', 'Non-protein Signaling']
    if any([pathway_type not in allowed_pathway_types for pathway_type in pathway_types]):
        raise ValueError(f"Invalid pathway type found in {pathway_types}. Please use one of the following: {allowed_pathway_types}.")
    
    senchat_results = adata.uns[senchat_output_key]
    masks = []

    if implausible_interactions:
        for sender, implausible_receivers in implausible_interactions.items():

            interactions_to_avoid = (senchat_results['sender'] == sender) & (senchat_results['receiver'].isin(implausible_receivers))

            if pathway_types:
                interactions_to_avoid &= (senchat_results['pathway_type'].isin(pathway_types))

            masks.append(~interactions_to_avoid)

        combined_mask = np.logical_and.reduce(masks) if masks else None

        constrained_results = senchat_results[combined_mask] if combined_mask is not None else senchat_results
    else: # Should be just paathway_types
        constrained_results = senchat_results[senchat_results['pathway_type'].isin(pathway_types)]
        
    adata.uns[f'{senchat_output_key}_constrained'] = constrained_results

def subset_for_communication(adata: AnnData,
                            senchat_output_key: str,
                            pval_threshold: float = None,
                            subset_sasp: bool = True,
                            layer: str = 'counts',
                            model: Literal['mouse', 'human'] = 'human') -> AnnData:


    data_dir = os.path.join(os.path.dirname(__file__), '..', 'data')

    if model not in ['human', 'mouse']:
        raise ValueError("Model must be either 'mouse' or 'human'. Working on expanding to more organisms.")
    
    data_path = os.path.join(data_dir, f'{model}_snc_genes.csv.gz')

    snc_genes = pd.read_csv(data_path, index_col=0)
    sasp_genes = snc_genes[snc_genes['SASP']].index.tolist()

    if senchat_output_key not in adata.uns:
        raise ValueError(f"Key '{senchat_output_key}' not found in adata.uns.")
    if layer not in adata.layers:
        raise ValueError(f"Layer '{layer}' not found in adata.layers.")
    
    senchat_output = adata.uns[senchat_output_key]

    if pval_threshold is not None:
        senchat_output = senchat_output[senchat_output['pval'] < pval_threshold]

    if subset_sasp:

        interactions_with_sasp_genes = senchat_output[senchat_output['ligand'].str.split('+').apply(lambda genes: any(gene in sasp_genes for gene in genes)) |
                                                 senchat_output['receptor'].str.split('+').apply(lambda genes: any(gene in sasp_genes for gene in genes))]

        senchat_ligs = list({gene for rec in interactions_with_sasp_genes['ligand'].unique() for gene in rec.split('+')})
        senchat_recs = list({gene for rec in interactions_with_sasp_genes['receptor'].unique() for gene in rec.split('+')})

    else:

        senchat_ligs = list({gene for rec in senchat_output['ligand'].unique() for gene in rec.split('+')})
        senchat_recs = list({gene for rec in senchat_output['receptor'].unique() for gene in rec.split('+')})

    ccc_genes = list(set(senchat_ligs + senchat_recs))

    adata_ccc = adata[:, ccc_genes]
    adata_ccc.X = adata_ccc.layers[layer]

    adata_ccc = adata_ccc[:, adata_ccc.X.sum(0) != 0]  
    adata_ccc = adata_ccc[adata_ccc.X.sum(1) != 0]

    print(f'Removed {adata.n_obs - adata_ccc.n_obs} cells with zero SASP counts.')
    print(f'Removed {len(ccc_genes)- adata_ccc.n_vars} genes from the initial SASP gene set.')
    print(f'SASP genes in the final dataset: {", ".join(ccc_genes)}')

    return adata_ccc

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

    print(top_genes.head(5))
    
    adata_cnmf.obsm['X_cNMF'] = usage.values
    adata_cnmf.varm['cNMF_spectra_scores'] = spectra_scores.values
    adata_cnmf.varm['cNMF_spectra_tpm'] = spectra_tpm.values
    # adata_cnmf.uns['cNMF_top_genes'] = top_genes

# Analysis (ligand enrichment), 
def calculate_ligand_enrichment(adata: AnnData,
                                spectra_tpm_key: str = 'cNMF_spectra_tpm') -> pd.DataFrame:

    if spectra_tpm_key not in adata.varm:
        raise ValueError(f"Key '{spectra_tpm_key}' not found in adata.varm.")
    
    spectra_tpm = adata.varm[spectra_tpm_key]

#  Calculate effect sizes (log-fold change, cohen's D)
# 
