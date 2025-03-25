import numpy as np
from anndata import AnnData
import pandas as pd
from itertools import product
import os

from collections.abc import Iterable as IterableClass
from typing import Dict, Literal, Optional, Sequence, Union, Callable, List
import contextlib
import joblib
from tqdm import tqdm

@contextlib.contextmanager
def tqdm_joblib(tqdm_object):
    """Context manager to patch joblib to report into tqdm progress bar given as argument."""
    class TqdmBatchCompletionCallback(joblib.parallel.BatchCompletionCallBack):
        def __call__(self, *args, **kwargs):
            tqdm_object.update(n=1)
            return super().__call__(*args, **kwargs)

    # Save original
    old_batch_callback = joblib.parallel.BatchCompletionCallBack
    joblib.parallel.BatchCompletionCallBack = TqdmBatchCompletionCallback
    try:
        yield tqdm_object
    finally:
        joblib.parallel.BatchCompletionCallBack = old_batch_callback
        tqdm_object.close()
from joblib import Parallel, delayed

def calculate_interaction_scores(adata: AnnData,
                                 sender_label: str,
                                 receiver_label: Optional[str] = None,
                                 sender_groups: Optional[List[str]] = None,
                                 receiver_groups: Optional[List[str]] = None,
                                 output_key: str = 'SenD_output',
                                 model: str = 'mouse',
                                 use_highly_variable: bool = False,
                                 highly_variable_key: Optional[str] = None,
                                 stringency: Optional[str] = 'neither',
                                 min_proportion: float = 0.1,
                                 test_permutation: bool = False,
                                 n_perms: int = 100,
                                 return_df: bool = False,
                                 n_jobs: int = 1) -> Optional[pd.DataFrame]:
    """ Calculate all interaction scores with respect to a specified lbel.

    Parameters
    ----------
    adata: AnnData
        Annotated data matrix containing expression data.
    group_label: str
        Key in `adata.obs` containing the group label.
    model: str
        The model organism. Options are 'mouse' or 'human'.
    use_highly_variable: bool
        Whether to use highly variable genes. If True, then you need to have previously determined which genes are highly variable, i.e.,
        the information shoudl be contained in adata.var[highly_variable_key].
    highly_variable_key: Optional[str]
        Key in adata.var containing the highly variable genes. If None, then the highly variable genes are not used.
    min_proportion: float
        Minimum proportion of cells expressing the ligand or receptor to be considered as a potential sender or receiver.
    test_permutation: bool
        Whether to calculate an associated p-value for the interaction scores using permutation testing.
    n_perms: int    
        Number of permutations to use for permutation testing
        
    Returns
    -------
    pd.DataFrame
        DataFrame containing the interaction scores across all ligand-receptor interaction pairs and potentially interacting groups
    """

    data_dir = os.path.join(os.path.dirname(__file__), '..', 'data')

    if model not in ['human', 'mouse']:
        raise ValueError("Model must be either 'mouse' or 'human'. Working on expanding to more organisms.")
    
    if use_highly_variable and stringency not in ['both', 'ligand', 'receptor', 'neither']:
            raise ValueError("stringency must be either 'both', 'ligand', 'receptor', or 'neither'.")
        
    data_path = os.path.join(data_dir, f'{model}_lr_pairs_and_tfs.csv.gz')
    all_lr_interactions = pd.read_csv(data_path, index_col=0)

    if use_highly_variable:
        if highly_variable_key not in adata.var:
            raise ValueError(f"Key {highly_variable_key} not found in adata.var.") 
          
        hvgs = adata.var_names[adata.var[highly_variable_key]]
        hvgs_set = set(hvgs) # Use set as it's faster formembership checking
    else:
        hvgs_set = None

    vars_set = set(adata.var_names) # Set is faster for checking

    if receiver_label is None:
        receiver_label = sender_label

    lr_pairs: Dict[tuple, dict] = {}
    for _, row in all_lr_interactions.iterrows():
        lig_subunits = row['ligand.symbol'].split(', ')
        rec_subunits = row['receptor.symbol'].split(', ')
        pathway_name = row['pathway_name']

        lig_subunits_set = set(lig_subunits)
        rec_subunits_set = set(rec_subunits)

        lig_tfs = []
        rec_tfs = []

        if pd.notna(row.get('Ligand-TF-combined')):
            lig_tfs = row['Ligand-TF-combined'].split('_')

        if pd.notna(row.get('Receptor-TF-combined')):
            rec_tfs = row['Receptor-TF-combined'].split('_')

        downstream_tfs_for_pair = list(set(lig_tfs + rec_tfs))

        # Filter out any TFs not in the dataset
        downstream_tfs_for_pair = [tf for tf in downstream_tfs_for_pair if tf in vars_set]

        # Only proceed if all ligands/receptors are in the data
        if not (lig_subunits_set.issubset(vars_set)\
            and rec_subunits_set.issubset(vars_set)):
            continue

        use_interaction = not use_highly_variable # If we're not using HVGs, then we should always the interaction
        if use_highly_variable:

            downstream_tfs_for_pair = [tf for tf in downstream_tfs_for_pair if tf in hvgs_set]

            # If we require all ligand sub-units and all receptor sub-units to be highly variable
            lig_in_hvg = lig_subunits_set.issubset(hvgs_set)
            rec_in_hvg = rec_subunits_set.issubset(hvgs_set)

            # If we require that all sub-units are highly variable
            if stringency == 'both':
                use_interaction = lig_in_hvg and rec_in_hvg 
                
            # If we only require that all ligand sub-units are highly variable
            elif stringency == 'ligand':
                use_interaction = lig_in_hvg and bool(lig_subunits_set & hvgs_set)
                
            # If we only require that all receptor sub-units are highly variable
            elif stringency == 'receptor':
                use_interaction = rec_in_hvg and bool(rec_subunits_set & hvgs_set)

            # If we only require that at least one ligand sub-unit and at least one-receptor subnit are in HVGs
            else:
                use_interaction = True

        if use_interaction:

            ligand = '+'.join(lig_subunits)
            receptor = '+'.join(rec_subunits)

            lr_pairs[(ligand, receptor)] = {'tf': downstream_tfs_for_pair, 'pathway_name': pathway_name}

    # Parallelise ligand-receptor score inference
    def _lr_score_for_one_pair(lig_str, rec_str, downstream_info):
        lig_list = lig_str.split('+')
        rec_list = rec_str.split('+')

        downstream_tfs = downstream_info['tf']
        pathway = downstream_info['pathway_name']

        return calculate_interaction_scores_for_lr_pair(adata,
                                                        ligand=lig_list,
                                                        receptor=rec_list,
                                                        sender_label=sender_label,
                                                        receiver_label=receiver_label,
                                                        sender_groups=sender_groups,
                                                        receiver_groups=receiver_groups,
                                                        downstream_tf=downstream_tfs or None,
                                                        pathway=pathway,
                                                        min_proportion=min_proportion,
                                                        test_permutation=test_permutation,
                                                        n_perms=n_perms)

    tasks = list(lr_pairs.items())  # each item is ((lig_str, rec_str), info)

    with tqdm_joblib(tqdm(desc='Computing ligand-receptor interaction scores', total=len(tasks))) as progress_bar:
        results = Parallel(n_jobs=n_jobs)(
            delayed(_lr_score_for_one_pair)(lig_str, rec_str, downstream_info)
            for (lig_str, rec_str), downstream_info in tasks
        )

    interaction_scores = pd.concat(results, ignore_index=True)

    stored_output_key = f'{output_key}_{sender_label}_{receiver_label}'

    if return_df:
        return interaction_scores
    else:
        adata.uns[stored_output_key] = interaction_scores
        print(f'Output stored in adata.uns[{stored_output_key}]')

def calculate_interaction_scores_for_lr_pair(adata: AnnData,
                                               ligand: Optional[Union[str, List[str]]],
                                               receptor: Optional[Union[str, List[str]]],
                                               sender_label: str,
                                               receiver_label: Optional[str] = None,
                                               sender_groups: Optional[List[str]] = None,
                                               receiver_groups: Optional[List[str]] = None,
                                               downstream_tf: Optional[Union[str, List[str]]] = None,
                                               pathway: Optional[str] = None,
                                               min_proportion: float =  0.1,
                                               test_permutation: bool = False,
                                               n_perms: int = 100
                                               ) -> pd.DataFrame:
    """ Calculate all interaction scores for a specified ligand-receptor pair 

    Parameters
    ----------
    adata: AnnData
        Annotated data matrix containing expression data.
    sender_label: str
        Key in `adata.obs` containing the group label.
    group1: str
        Name of the first group.
    group2: str
        Name of the second group.
    ligand_receptor_pairs: pd.DataFrame
        DataFrame containing ligand-receptor pairs.
    interaction_score_method: str
        Method to calculate the interaction score. Options are "geometric_mean" or "arithmetic_mean".
    tf_expression: Optional[Union[str, np.ndarray]]
        Key in `adata.obs` containing the expression of a downstream transcription factor. If None, the interaction score is calculated without considering the downstream transcription factor.

    Returns
    -------
    pd.DataFrame
        DataFrame containing the interaction scores for a specified ligand-receptor pair
        between all potentially interacting groups.
    """

    if sender_groups is None:
        sender_groups = adata.obs[sender_label].unique()

    if receiver_groups is None:
        
        if receiver_label is None:
            receiver_label = sender_label
            
        receiver_groups = adata.obs[receiver_label].unique()

    # Pre-compute sender and receiver group indices to save time
    sender_indices_map = {
        group: np.where(adata.obs[sender_label] == group)[0]
        for group in sender_groups
    }
    receiver_indices_map = {
        group: np.where(adata.obs[receiver_label] == group)[0]
        for group in receiver_groups
    }

    ligand_expression = adata[:, ligand].X.toarray()
    receptor_expression = adata[:, receptor].X.toarray()

    # Pre-compute the cell-wise ligand and receptor geometric means to speed up permutations
    with np.errstate(divide='ignore', invalid='ignore'):
        ligand_geom_mean = np.exp(np.log(ligand_expression).mean(1))
        receptor_geom_mean = np.exp(np.log(receptor_expression).mean(1))

    tf_expression = None
    if downstream_tf is not None:
        tf_expression = adata[:, downstream_tf].X.mean(1).A1

    # We define a potential sender as a group as any group with at least min_proportion of its cells expressing the ligand.
    # We define a potential receiver similarly
    potential_sender_receiver_pairs = []
    for sender, receiver in product(sender_groups, receiver_groups):

        sender_indices = sender_indices_map[sender]
        receiver_indices =  receiver_indices_map[receiver]
 
        proportion_expressing_ligand = (ligand_expression[sender_indices] > 0).all(axis=1).mean() # The geometric mean being > 0 is equivalent to all of them being > 0

        proportion_expressing_receptor = (receptor_expression[receiver_indices] > 0).all(axis=1).mean()
        
        if ( (proportion_expressing_ligand > min_proportion)\
            &(proportion_expressing_receptor > min_proportion) ):

            potential_sender_receiver_pairs.append((sender, receiver))

    interaction_score_results = []
    # Calculate the scores per sender-receiver pair
    for (sender, receiver) in potential_sender_receiver_pairs:

        sender_indices = sender_indices_map[sender]
        receiver_indices = receiver_indices_map[receiver]

        ligand_expression_in_sender = ligand_geom_mean[sender_indices]
        receptor_expression_in_receiver = receptor_geom_mean[receiver_indices]

        tf_expression_in_receiver = tf_expression[receiver_indices] if tf_expression is not None else None

        ligand_mean = ligand_expression_in_sender.mean()
        receptor_mean = receptor_expression_in_receiver.mean()
        interaction_score = ligand_mean * receptor_mean

        # If we specify downstream TF expression, we take the arithmetic mean, as we only need one to be "signficantly" activated
        if tf_expression is not None:
                
            tf_mean = tf_expression_in_receiver.mean()
            interaction_score = ligand_mean * (receptor_mean * tf_mean)**(0.5)

        pval = np.nan
        if test_permutation:
            
            n_senders = len(sender_indices)
            n_receivers = len(receiver_indices)
            all_indices = np.arange(adata.n_obs)

            interaction_scores_null = np.zeros(n_perms)
            for perm in range(n_perms):

                sender_indices_perm = np.random.choice(all_indices, n_senders, replace=False) # Sample sender group
                receiver_indices_perm = np.random.choice(all_indices, n_receivers, replace=False) # Sample receiver group

                ligand_expression_in_sender_perm = ligand_geom_mean[sender_indices_perm]
                receptor_expression_in_receiver_perm = receptor_geom_mean[receiver_indices_perm]

                tf_expression_in_receiver_perm = tf_expression[receiver_indices_perm] if tf_expression is not None else None

                ligand_mean_perm = ligand_expression_in_sender_perm.mean()
                receptor_mean_perm = receptor_expression_in_receiver_perm.mean()
                interaction_score_perm =  ligand_mean_perm * receptor_mean_perm

                # If we specify downstream TF expression, we take the arithmetic mean, as we only need one to be "signficantly" activated
                if tf_expression is not None:
                        
                    tf_mean_perm = tf_expression_in_receiver_perm.mean()
                    interaction_score_perm = ligand_mean_perm * (receptor_mean_perm * tf_mean_perm)**(0.5)

                interaction_scores_null[perm] = interaction_score_perm

            pval = (interaction_score < interaction_scores_null).sum() / n_perms
            
        # Store
        interaction_score_results.append({
            'sender': sender,
            'receiver': receiver,
            'interaction_score': interaction_score,
            'ligand': '+'.join(sorted(ligand)) if isinstance(ligand, list) else ligand,
            'receptor': '+'.join(sorted(receptor)) if isinstance(receptor, list) else receptor,
            'interaction': f"{'+'.join(sorted(ligand)) if isinstance(ligand, list) else ligand}"
                           f" - "
                           f"{'+'.join(sorted(receptor)) if isinstance(receptor, list) else receptor}",
            'downstream_tfs': (
                ', '.join(sorted(downstream_tf))
                if downstream_tf is not None and isinstance(downstream_tf, list) else ''
            ),
            'pathway': pathway,
            'pval': pval
        })

    # Convert to DataFrame
    return pd.DataFrame(interaction_score_results)
