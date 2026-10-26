

from typing import Tuple, Optional
import monod
from itertools import product
import anndata as ad 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from monod.analysis import DE_parameters
from monod.extract_data import extract_data
from monod.inference import perform_inference, reject_genes
from monod.cme_toolbox import CMEModel

# def _get_log_fold_changes(sr1: inference.SearchResults,
#                           sr2: inference.SearchResults,
#                           sd1: extract_data.SearchData,
#                           sd2: extract_data.SearchData,
#                           dname1: str,
#                           dname2: str,
#                           gene_names: list[str],
#                           gf_rej: bool = False,
#                           param_lfc: float = 2.0,
#                           mean_lfc: float = 1.0,
#                           pval_thr: float = 0.05,
#                           outlier_de: bool = True,
#                           single_nuc: bool = False,
#                           correct_off: bool = False):
#     '''
#     Utilize different metrics to find fold-changes (FCs) between clusters in different SRs

#     sr1: SearchResults object 1 (single object)
#     sr2: SearchResults object 2
#     sd1: SearchData object for sr1
#     sd2: SearchData object for sr2
#     gene_names: list of gene names to compare between SRs
#     gf_rej: whether to use boolean list of rejected genes from both SRs
#     param_lfc: FC threshold value (to call DE-theta genes)
#     mean_lfc: Mean S expression threshold value, for genes to consider
#     outlier_de: Use iterative outlier calling procedure to assign DE-theta genes (see Monod https://github.com/pachterlab/monod_examples/blob/main/Monod_demo.ipynb)
#     single_nuc: is this single_nuclear RNA data
#     correct_off: boolean to correct parameter offset with ODR
#     '''

#     param_lfcs = pd.DataFrame()
#     fcs,types,which_pair,highFC,spliceFC,g_names,out_de = ([] for i in range(7))

#     sr1 = sr1
#     sr2 = sr2

#     gfilt1 = [list(sr1.gene_names).index(i) for i in gene_names]
#     gfilt2 = [list(sr2.gene_names).index(i) for i in gene_names]

#     par_vals1 = np.copy(sr1.param_estimates[:,gfilt1,:])
#     par_vals2 = np.copy(sr2.param_estimates[:,gfilt2,:])

#     parnames = ('b','beta','gamma')

#     if correct_off:
#         param_names = sr1.model.get_log_name_str()
#         offsets = []
#         for k in range(len(param_names)):
#             m1 = par_vals1[0,:,k]
#             m2 = par_vals2[0,:,k]
#             offset = analysis.diffexp_fpi(m1,m2,param_names[k],viz=False, pval_thr=pval_thr)[1]
#             par_vals2[0,:,k] -= offset

#         fc_par = (par_vals1-par_vals2)/np.log10(2)
#     else:
#         fc_par = (par_vals1-par_vals2)/np.log10(2)  #Get FCs between cluster params

#     print('fc_par.shape: ',fc_par.shape)
#     if single_nuc:
#         fc_s_par = np.log2((sd1.layers[0][gfilt1,:].mean(1) + 1e-12)/(sd2.layers[0][gfilt2,:].mean(1) + 1e-12)) #Get unspliced FCs
#     else:
#         fc_s_par = np.log2((sd1.layers[1][gfilt1,:].mean(1) + 1e-12)/(sd2.layers[1][gfilt2,:].mean(1) + 1e-12)) #Get spliced FCs

#     print('fc_s_par.shape: ',fc_s_par.shape)

#     #Use outlier calling method to find DE genes
#     if outlier_de:
#         if len(sr1.gene_names) != len(sr2.gene_names):
#             print('Not running outlier DE. SRs need to have the same gene dimensions.')
#             par_bool_de = np.zeros((len(gene_names),len(parnames)))
#         else:
#             dr_analysis = monod.analysis.diffexp_pars(sr1,sr2,viz=True,modeltype='id',use_sigma=True, pval_thr=pval_thr)
#             par_bool_de = dr_analysis[1].T

#   #-----is parameter FC significant -----
#     if gf_rej is False:
#         gf_rej = [True]*(len(gfilt1))
#     else:
#         gf_rej = ~(sr1.rejected_genes[gfilt1]) & ~(sr2.rejected_genes[gfilt2])


#     for n in range(len(parnames)):
#         #Boolean for if large param FC and not rejected gene (with minimum expression)
#         if single_nuc:
#             gf_highnoise = (np.abs(fc_par[0,:,n])>param_lfc)  \
#                 & ((sd1.layers[0][gfilt1,:].mean(1)>mean_lfc) | (sd2.layers[0][gfilt2,:].mean(1)>mean_lfc)) \
#                 & gf_rej
#         else:
#             gf_highnoise = (np.abs(fc_par[0,:,n])>param_lfc)  \
#                 & ((sd1.layers[1][gfilt1,:].mean(1)>mean_lfc) | (sd2.layers[1][gfilt2,:].mean(1)>mean_lfc)) \
#                 & gf_rej

#         #Boolean for FC (above) but no FC detected at S-level
#         gf_highnoise_meanS = gf_highnoise & (np.abs(fc_s_par)<mean_lfc) & gf_rej

#         #Boolean for FC (above)
#         gf_onlyhigh = gf_highnoise & gf_rej

#         #For dataframe
#         fcs += list(fc_par[0,gf_rej,n])
#         g_names += list(gene_names[gf_rej])
#         which_pair += [[dname1, dname2]]*np.sum(gf_rej)
#         highFC += list(gf_onlyhigh[gf_rej])
#         spliceFC += list(gf_highnoise_meanS[gf_rej])
#         types += [parnames[n]]*np.sum(gf_rej)
#         if outlier_de:
#             out_de += list(par_bool_de[gf_rej,n])

#     if outlier_de:
#         param_lfcs['deTheta_outlier'] = out_de

#     param_lfcs['log2FC'] = fcs
#     param_lfcs['gene'] = g_names
#     param_lfcs['cluster_pair'] = which_pair
#     param_lfcs['deTheta_FC'] = highFC
#     param_lfcs['deTheta_noDeMuS'] = spliceFC
#     param_lfcs['param'] = types

#     return param_lfcs

# def _get_DE_results(param_lfcs: pd.DataFrame) -> pd.DataFrame:

#     # Pre-filtering by parameter
#     df_b = param_lfcs[param_lfcs.param == 'b'].reset_index(drop=True)
#     df_beta = param_lfcs[param_lfcs.param == 'beta'].reset_index(drop=True)
#     df_gamma = param_lfcs[param_lfcs.param == 'gamma'].reset_index(drop=True)

#     # Extracting vectors
#     bs = df_b.log2FC
#     betas = df_beta.log2FC
#     gammas = df_gamma.log2FC
#     names = df_b.gene

#     # Combining significance masks
#     sig_b = df_b.deTheta_FC | df_b.deTheta_noDeMuS
#     sig_beta = df_beta.deTheta_FC | df_beta.deTheta_noDeMuS
#     sig_gamma = df_gamma.deTheta_FC | df_gamma.deTheta_noDeMuS

#     highFCs = df_b.deTheta_FC.values | df_beta.deTheta_FC.values | df_gamma.deTheta_FC.values
#     noSpliceFCs = df_b.deTheta_noDeMuS.values | df_beta.deTheta_noDeMuS.values | df_gamma.deTheta_noDeMuS.values

#     dom = []
#     dom_type = []

#     # Prepare split cluster pairs
#     pairs = df_b.cluster_pair.apply(lambda x: x if isinstance(x, list) else x.split(','))
#     pairs = np.vstack(pairs.values)
#     # pairs = df_b.cluster_pair.str.split(',').tolist()
#     # pairs = np.array(pairs)

#     # Precompute necessary comparisons
#     abs_b = np.abs(bs)
#     abs_beta = np.abs(betas)
#     abs_gamma = np.abs(gammas)

#     # Vectorised way of considering significance conditions
#     cond1 = (-betas > 0) & sig_beta & (-gammas > 0) & sig_gamma # Possibly significant in k (sig. change in both beta and k)
#     cond2 = (betas > 0) & sig_beta & (gammas > 0) & sig_gamma # Possibly significant in k (sig. change in both beta and k)
#     cond3 = sig_b & (abs_b > abs_beta) & (abs_b > abs_gamma) # Sig in b
#     cond4 = sig_b & (abs_b < abs_beta) & (abs_b >= abs_gamma) # Sig in b
#     cond5 = sig_b & (abs_b < abs_gamma) & (abs_b >= abs_beta) # Sig in b
#     cond6 = (betas > 0) & sig_beta
#     cond7 = (gammas > 0) & sig_gamma # Sig in gamma
#     cond8 = (betas < 0) & sig_beta
#     cond9 = (gammas < 0) & sig_gamma # Sig in gamma

#     # Initialize outputs
#     dom = np.full(len(bs), 'None', dtype='<U50')
#     dom_type = np.full(len(bs), 'None', dtype='<U50')

#     # Apply conditions vectorized
#     mask = cond1 & (abs_b < abs_beta)
#     dom[mask] = pairs[mask, 0]
#     dom_type[mask] = 'frequency'

#     mask = cond1 & (bs > 0) & sig_b
#     dom[mask] = pairs[mask, 0]
#     dom_type[mask] = 'frequency'

#     mask = cond1 & (bs < 0) & sig_b
#     dom[mask] = pairs[mask, 1]
#     dom_type[mask] = 'frequency'

#     mask = cond2 & (abs_b < abs_beta)
#     dom[mask] = pairs[mask, 1]
#     dom_type[mask] = 'frequency'

#     mask = cond2 & (bs < 0) & sig_b
#     dom[mask] = pairs[mask, 1]
#     dom_type[mask] = 'frequency'

#     mask = cond3 & (bs > 0)
#     dom[mask] = pairs[mask, 0]
#     dom_type[mask] = 'synthesis'

#     mask = cond3 & (bs < 0)
#     dom[mask] = pairs[mask, 1]
#     dom_type[mask] = 'synthesis'

#     mask = cond4 & (betas > 0)
#     dom[mask] = pairs[mask, 0]
#     dom_type[mask] = 'synthesis'

#     mask = cond4 & (betas < 0)
#     dom[mask] = pairs[mask, 1]
#     dom_type[mask] = 'synthesis'

#     mask = cond5 & (gammas < 0)
#     dom[mask] = pairs[mask, 0]
#     dom_type[mask] = 'synthesis'

#     mask = cond5 & (gammas > 0)
#     dom[mask] = pairs[mask, 1]
#     dom_type[mask] = 'synthesis'

#     # Beta/gamma individually if none of the above
#     mask = cond6 & ~cond1 & ~cond2 & ~sig_b
#     dom[mask] = pairs[mask, 0]
#     dom_type[mask] = 'splicing'

#     mask = cond7 & ~cond1 & ~cond2 & ~sig_b & ~cond6
#     dom[mask] = pairs[mask, 1]
#     dom_type[mask] = 'degradation'

#     mask = cond8 & ~cond1 & ~cond2 & ~sig_b & ~cond6 & ~cond7
#     dom[mask] = pairs[mask, 1]
#     dom_type[mask] = 'splicing'

#     mask = cond9 & ~cond1 & ~cond2 & ~sig_b & ~cond6 & ~cond7 & ~cond8
#     dom[mask] = pairs[mask, 0]
#     dom_type[mask] = 'degradation'

#     # Approximate the log fold change in k using the beta and gamma log fold changes
#     ks = -0.5 * (betas.values + gammas.values)

#     markers_df = pd.DataFrame(data={'gene': names.values,
#                               'synthesis': bs.values,
#                                 'splicing': -betas.values,
#                                 'degradation': -gammas.values,
#                                 'frequency': ks,
#                                 'DE': highFCs,
#                                 'DE-reg': noSpliceFCs,
#                                 'marker': dom,
#                                 'marker_type': dom_type,
#                                 })

#     markers_df['marker'] = pd.Categorical(markers_df['marker'], categories=np.unique(pairs), ordered=False)
#     markers_df['marker_type'] = pd.Categorical(markers_df['marker_type'],
#                                    categories=['None', 'synthesis', 'frequency', 'splicing', 'degradation'],
#                                    ordered=False)

#     return markers_df

def prepare_files_for_inference(loom_directory: str,
                                dataset_names: list[str],
                                transcriptome_filepath: str
                                ) -> Tuple[str, list[str]]: 

    loom_filepaths = [f'{loom_directory}/{dname}.loom' for dname in dataset_names]

    adata = ad.read_loom(loom_filepaths[0])
    genes_to_use = list(adata.var['gene_name'].values)
    n_genes = len(genes_to_use)

    return monod.preprocess.construct_batch(dataset_filepaths=loom_filepaths,
                                            transcriptome_filepath=transcriptome_filepath,
                                            dataset_names=dataset_names,
                                            attribute_names=[('unspliced','spliced'), 'gene_name', 'barcode'],
                                            meta='',
                                            creator='SenID',
                                            datestring='',
                                            exp_filter_threshold=None,
                                            n_genes=n_genes,
                                            genes_to_fit=genes_to_use)

def run_monod_inference(adata: ad.AnnData,
                        transcriptome_filepath: str,
                        dataset_str: str,
                        fitmodel: CMEModel = CMEModel('Bursty','Poisson'),
                        genes_to_fit: Optional[list] = None,
                        **kwargs
                        ) -> ad.AnnData:

   # Determine n_genes based on the filter.
    n_genes = len(genes_to_fit) if genes_to_fit is not None else adata.n_vars
                            
    return perform_inference(adata,
                            transcriptome_filepath=transcriptome_filepath,
                            fitmodel=fitmodel,
                            dataset_string = dataset_str,
                            genes_to_fit = genes_to_fit,
                            **kwargs)

def process_monod_fits(adata: ad.AnnData,
                       **kwargs
                       ) -> ad.AnnData:

    return reject_genes(adata, **kwargs)

def get_monod_de_genes(adata1: ad.AnnData,
                       adata2: ad.AnnData,
                       **kwargs) -> None:

    _, _, _, _ = DE_parameters(adata1,
                                adata2,
                                **kwargs)
                                                    
