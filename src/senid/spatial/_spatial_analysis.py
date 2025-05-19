
from typing import Literal, Dict
import tensorflow as tf  # For building and training the machine learning model
import pandas as pd  # For data manipulation and analysis
from anndata import AnnData
from scanpy import pp  # For preprocessing single-cell data
import scanpy as sc
import random
import os
from os import path  # For pathname manipulations
from mNSF import MoranI # For calculating spatial dependeicy of each factor within each sample, using Moran's I
from itertools import product
import numpy as np
import warnings

# Ignore all DeprecationWarnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

import mNSF
from mNSF import process_multiSample
from mNSF.NSF import preprocess, misc, visualize
from mNSF import training_multiSample

# nchunk = 1
# nsample = 4

# # Set up paths for outputs
# mpth = path.join("models")
# misc.mkdir_p(mpth)
# pp = path.join(mpth, "pp", str(2))
# misc.mkdir_p(pp)

# list_D = []
# list_X = []

# data_directory = '/Users/axelalmet/Documents/CellCellCommunicationModelling/Data/Park2024/'
# samples = ['P1', 'P2', 'P3', 'P4']

# for ksample in range(nsample):
#     adata_sample = sc.read_h5ad(data_directory + 'park24_' + samples[ksample] + '_sasp.h5ad')
#     Y = pd.DataFrame(adata_sample.layers['counts'].toarray())
#     X = pd.DataFrame(adata_sample.obsm['spatial'])
#     D = process_multiSample.get_D(X, Y)
#     list_D.append(D)
#     list_X.append(D["X"])

# # list_sampleID = process_multiSample.get_listSampleID(list_D)
# list_Dtrain = process_multiSample.get_listDtrain(list_D)

# # get goodness-of-fit poisson deviance for each sample for each factor
# Lvals = range(2, 21)
# deviance_values = np.zeros(len(Lvals))

# for i, L in enumerate(Lvals):

#     list_fit = process_multiSample.ini_multiSample(list_D, L, "nb")
#     vec_dev = 0
#     for ksample in range(0,nsample):
#         dev_mnsf = visualize.gof(list_fit[ksample],list_D[ksample],Dval=None,S=10,plot=False)
#         vec_dev= vec_dev + dev_mnsf['tr']['mean']

#     deviance_values[i] = vec_dev
#     # print("L="+str(L))
#     # print("deviance="+str(vec_dev))
#     # print("")

# plt.plot(Lvals, deviance_values)
# plt.savefig(data_directory +  'park24_mnsf_poisson_deviance_values.pdf', bbox_inches='tight')
# plt.show()

# L = 10 # Gave best Poisson deviance fit, the others didn't really improve significantly 
	
# # Process data chunking
# list_D_chunked=list()
# list_X_chunked=list()
# list_chunk_mapping = list()
# for ksample in range(0,nsample):
#     adata_sample = sc.read_h5ad(data_directory + 'park24_' + samples[ksample] + '_sasp.h5ad')
#     Y = pd.DataFrame(adata_sample.layers['counts'].toarray())
#     X = pd.DataFrame(adata_sample.obsm['spatial'])
#     list_D_sampleTmp,list_X_sampleTmp, chunk_mapping = process_multiSample.get_chunked_data(X.iloc[:,:],Y.iloc[:,:],nchunk,method = "random") #choose method = "balanced_kmeans" for chunking the spots based on the spatial coordinates
#     list_D_chunked = list_D_chunked + list_D_sampleTmp
#     list_X_chunked = list_X_chunked + list_X_sampleTmp
#     list_chunk_mapping.append(chunk_mapping)

# print(list_chunk_mapping)
# # #  Extracts the training data from our processed data. This function prepares the data in the format required for model training.
# list_Dtrain = process_multiSample.get_listDtrain(list_D_chunked)
# list_sampleID=process_multiSample.get_listSampleID(list_D)

# # Set up induced points for each sample
# for ksample in range(nsample):
#     # Select 15% of spots as induced points
#     ninduced = round(list_D_chunked[ksample]['X'].shape[0] * 0.35)
#     rd_ = random.sample(range(list_D_chunked[ksample]['X'].shape[0]), ninduced)
#     list_D_chunked[ksample]["Z"] = list_D_chunked[ksample]['X'][rd_, :]

# list_fit = process_multiSample.ini_multiSample(list_D_chunked, L, "nb")

# list_fit = training_multiSample.train_model_mNSF(
#     list_fit,      # Initialized model
#     pp,            # Directory for preprocessing results
#     list_Dtrain,    # Chunked training data
#     list_D_chunked, # Full chunked dataset
#     num_epochs=500,  # Number of training iterations
#     nsample = nsample, # Number of samples
#     nchunk = nchunk, # Number of chunks
#     legacy=True, # Use legacy mode for training as working with M1/M2 mac
# )

# process_multiSample.save_object(list_fit, 'park24_sasp_list_fit_nb_' + str(nsample) + 'samples_szMean_L' + str(L) + '_fullData.pkl') 

# # Load object
# import pickle
# with open('park24_sasp_list_fit_nb_' + str(nsample) + 'samples_szMean_L' + str(L) + '_fullData.pkl', 'rb') as f:
#     list_fit = pickle.load(f)

# inpf12=process_multiSample.interpret_npf_v3(list_fit,list_X,S=100,lda_mode=False)


# W = inpf12["loadings"]
# #Wdf=pd.DataFrame(W*inpf12["totals1"
# Wdf=pd.DataFrame(W*inpf12["totalsW"][:,None],  columns=range(1,L+1))
# Wdf.to_csv(path.join('park24_sasp_loadings_spde_nb_' +  str(nsample) + 'samples_szMean_L' + str(L) + '_fullDataa.csv'))

# # ## save the factors
# Factors = inpf12["factors"][:,0:L]

# for k in range(0,nsample):
# 	indices=list_ s

# Fplot = misc.t2np(list_fit[0].sample_latent_GP_funcs(list_D_chunked[0]["X"], S=3, chol=False)).T

# # print(Fplot.shape)
# # print(list_X)
# # print(list_chunk_mapping)

# # print(sum(len(X) for X in list_X))
# Fplot = process_multiSample.reorder_chunked_results(Fplot,list_chunk_mapping,list_X)

# for i in range(L):
#     I, p_value = MoranI.calculate_morans_i(list_D[0]["X"], Fplot[:, i])
#     print(f"Factor {i+1} - Moran's I: {I:.4f}, p-value: {p_value:.4f}")

def subset_for_spatial_sasp_communication(adata: AnnData | Dict[str],
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
