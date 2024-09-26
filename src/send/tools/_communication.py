import numpy as np
from anndata import AnnData
import pandas as pd
from itertools import product

from collections.abc import Iterable as IterableClass
from typing import Dict, Literal, Optional, Sequence, Union, Callable, List

def calculate_interaction_scores(adata: AnnData,
                                 group_label: str,
                                 model: str = 'mouse',
                                 use_highly_variable: bool = False,
                                 highly_variable_key: Optional[str] = None,
                                 highly_variable_stringency: Optional[str] = 'stringent_both',
                                 min_proportion: float = 0.1,
                                 test_permutation: bool = False,
                                 n_perms: int = 100) -> pd.DataFrame:
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

    if model not in ['human', 'mouse']:
        raise ValueError("Model must be either 'mouse' or 'human'. Working on expanding to more organisms.")
    
    if use_highly_variable:
        if highly_variable_key is None:
            raise ValueError("If use_highly_variable is True, highly_variable_key must be specified.")
    
        if highly_variable_key not in adata.var.keys():
            raise ValueError(f"Key {highly_variable_key} not found in adata.var.keys().")   
    
        if highly_variable_stringency not in ['both', 'ligand', 'receptor', 'neither']:
            raise ValueError("highly_variable_stringency must be either 'both', 'ligand', 'receptor', or 'neither'.")
        
    all_lr_interactions = pd.read_csv(f'data/{model}_lr_pairs_and_tfs.csv', index_col=0)

    lr_pairs = {}
    for i, row in all_lr_interactions.iterrows():
        ligand_symbol = row['ligand.symbol'].split(', ')
        receptor_symbol = row['receptor.symbol'].split(', ')
        pathway_name = row['pathway_name']

        lig_tfs = []
        rec_tfs = []

        if 'Ligand-TF-combined' in row.dropna():
            lig_tfs = row['Ligand-TF-combined'].split('_')

        if 'Receptor-TF-combined' in row.dropna():
            rec_tfs = row['Receptor-TF-combined'].split('_')

        downstream_tfs_for_pair = list(set(lig_tfs + rec_tfs))

        # Filter out any TFs not in the dataset and not in HVGs
        downstream_tfs_for_pair = [tf for tf in downstream_tfs_for_pair if tf in adata.var_names]

        use_interaction = False
        if set(ligand_symbol).issubset(adata.var_names)\
            and set(receptor_symbol).issubset(adata.var_names):


            if use_highly_variable:

                downstream_tfs_for_pair = [tf for tf in downstream_tfs_for_pair if tf in adata.var_names[adata.var[highly_variable_key]]]

                # If we require all ligand sub-units and all receptor sub-units to be highly variable
                if highly_variable_stringency == 'both':
                    if set(ligand_symbol).issubset(set(adata.var_names[adata.var[highly_variable_key]]))\
                        and set(receptor_symbol).issubset(set(adata.var_names[adata.var[highly_variable_key]])):
                        use_interaction = True

                # If we only require that all ligand sub-units are highly variable
                elif highly_variable_stringency == 'ligand':
                    if set(ligand_symbol).issubset(set(adata.var_names[adata.var[highly_variable_key]]))\
                        and bool(set(receptor_symbol) & set(set(adata.var_names[adata.var[highly_variable_key]]))):
                        use_interaction = True
                    
                # If we only require that all receptor sub-units are highly variable
                elif highly_variable_stringency == 'receptor':
                    if bool(set(ligand_symbol) & set(set(adata.var_names[adata.var[highly_variable_key]])))\
                        and set(receptor_symbol).issubset(set(adata.var_names[adata.var[highly_variable_key]])):
                        use_interaction = True

                # If we only require that at least one ligand sub-unit and at least one-receptor subnit are in HVGs
                else:
                    use_interaction = True

            if use_interaction:

                ligand = row['ligand.symbol'].replace(', ', '+')
                receptor = row['receptor.symbol'].replace(', ', '+')

                lr_pairs[(ligand, receptor)] = {'tf': downstream_tfs_for_pair, 'pathway_name': pathway_name}

    lr_interaction_scores = []
    for pair in lr_pairs:

        ligand_unit, receptor_unit = pair
        ligand = ligand_unit.split('+')
        receptor = receptor_unit.split('+')

        downstream_tfs = lr_pairs[pair]['tf']
        pathway = lr_pairs[pair]['pathway_name']

        ligand_expression = adata[:, ligand].X.toarray()
        receptor_expression = adata[:, receptor].X.toarray()

        if len(downstream_tfs) != 0:

            interaction_score = calculate_interaction_scores_for_lr_pair(adata,
                                                                        group_label,
                                                                        ligand,
                                                                        receptor,
                                                                        downstream_tfs,
                                                                        pathway,
                                                                        min_proportion,
                                                                        test_permutation,
                                                                        n_perms)
            
        else:
            interaction_score = calculate_interaction_scores_for_lr_pair(adata,
                                                                        group_label,
                                                                        ligand,
                                                                        receptor,
                                                                        pathway,
                                                                        min_proportion,
                                                                        test_permutation,
                                                                        n_perms)

        lr_interaction_scores.append(interaction_score)

    interaction_scores = pd.concat(lr_interaction_scores)
    
    return interaction_scores

def calculate_interaction_scores_for_lr_pair(adata: AnnData,
                                               ligand: Optional[Union[str, List[str]]],
                                               receptor: Optional[Union[str, List[str]]],
                                               sender_label: str,
                                               receiver_label: Optional[str] = None,
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

    sender_groups = adata.obs[sender_label].unique()

    if receiver_label is None:
        receiver_groups = adata.obs[sender_label].unique()
    else:
        receiver_groups = adata.obs[receiver_label].unique()

    ligand_expression = adata[:, ligand].X.toarray()
    receptor_expression = adata[:, receptor].X.toarray()

    if downstream_tf is not None:
        tf_expression = adata[:, downstream_tf].X.toarray()

    # We define a potential sender as a group as any group with at least min_proportion of its cells expressing the ligand.
    # We define a potential receiver similarly
    potential_sender_receiver_pairs = []
    for sender, receiver in product(sender_groups, receiver_groups):

        sender_indices = np.where(adata.obs[sender_label] == sender)[0]
        receiver_indices =  np.where(adata.obs[receiver_label] == receiver)[0]

        proportion_expressing_ligand = (np.exp(np.log(ligand_expression[sender_indices]).sum(1)) > 0).sum()\
                                        / len(sender_indices)

        proportion_expressing_receptor = (np.exp(np.log(receptor_expression[receiver_indices]).sum(1)) > 0).sum()\
                                        / len(receiver_indices)
        
        if (proportion_expressing_ligand > min_proportion)\
            &(proportion_expressing_receptor > min_proportion)\
            &((sender, receiver) not in potential_sender_receiver_pairs):

            potential_sender_receiver_pairs.append((sender, receiver))

    interaction_scores = []
    pvals = []
    # Calculate the scores per sender-receiver pair
    for pair in potential_sender_receiver_pairs:

        sender, receiver = pair
        sender_indices = np.where(adata.obs[sender_label] == sender)[0]
        receiver_indices = np.where(adata.obs[receiver_label] == receiver)[0]

        ligand_expression_in_sender = ligand_expression[sender_indices]
        receptor_expression_in_receiver = receptor_expression[receiver_indices]

        if downstream_tf is not None:
            tf_expression_in_receiver = tf_expression_in_receiver[receiver_indices]

        interaction_score = calculate_interaction_score_for_lr_pair(ligand_expression_in_sender,
                                                                    receptor_expression_in_receiver,
                                                                    tf_expression_in_receiver)
        
        interaction_scores.append(interaction_score)

        if test_permutation:
            
            n_senders = len(sender_indices)
            n_receivers = len(receiver_indices)

            interaction_scores_null = np.zeros(n_perms)

            possible_sender_indices = np.arange(adata.n_obs)
            possible_receiver_indices = np.arange(adata.n_obs)


            for perm in range(n_perms):

                # Sample sender group
                sender_indices_perm = np.random.choice(possible_sender_indices, n_senders, replace=False)

                # Sample receiver group
                receiver_indices_perm = np.random.choice(possible_receiver_indices, n_receivers, replace=False)

                ligand_expression_in_sender_perm = ligand_expression[sender_indices_perm]
                receptor_expression_in_receiver_perm = receptor_expression[receiver_indices_perm]

                if downstream_tf is not None:
                    tf_expression_in_receiver_perm = tf_expression_in_receiver[receiver_indices_perm]

                interaction_score_perm = calculate_interaction_score_for_lr_pair(ligand_expression_in_sender_perm,
                                                                                receptor_expression_in_receiver_perm,
                                                                                tf_expression_in_receiver_perm)
                

                interaction_scores_null[perm] = interaction_score_perm

            pval_for_pair = (interaction_score < interaction_scores_null).sum() / n_perms
            pvals.append(pval_for_pair)

    potential_senders = [pair[0] for pair in potential_sender_receiver_pairs]
    potential_receivers = [pair[1] for pair in potential_sender_receiver_pairs] 

    interactions_for_lr_pair = pd.DataFrame({'sender': potential_senders,
                                             'receiver': potential_receivers,
                                             'interaction_score': interaction_scores})
    
    # Specify the information about the ligand-receptor interaction
    interactions_for_lr_pair['ligand'] = '+'.join(ligand)
    interactions_for_lr_pair['receptor'] = '+'.join(receptor)

    interactions_for_lr_pair['interaction'] = '+'.join(ligand) + '-' + '+'.join(receptor)
    interactions_for_lr_pair['downstream_tfs'] = ''

    if  downstream_tf is not None:
        interactions_for_lr_pair['downstream_tfs'] = '+'.join(downstream_tf)

    interactions_for_lr_pair['pathway'] = pathway  

    if test_permutation:
        interactions_for_lr_pair['pval'] = pvals                                      

    return interactions_for_lr_pair

def calculate_interaction_score_for_lr_pair(ligand_expression: np.ndarray,
                                receptor_expression: np.ndarray,
                                tf_expression: np.ndarray = None,
                                ) -> float:
    """ Calculate the interaction score between a ligand and receptor pair. Formula is based on a simple multiplicative formula,
    I = L*R*TF, where, in the case of ligand or receptor multi-units, L and R are taken as geometric means. If downstream TFs are specified,
    TF is taken as the arithmetic mean. 
    
    Parameters
    ----------
    ligand_expression: np.ndarray
        The ligand expression. May be 1D or ND to account for ligand multi-units.
    receptor_expression: np.ndarray
        The receptor expression. May be 1D or ND to account for ligand multi-units.
  
    Returns
    -------
    np.ndarray
        The interaction score between the two groups of cells.
    """

    # This expression accounts for the fact that we take the geometric mean when calculating ligand and receptor expression
    ligand_mean = np.exp(np.log(ligand_expression).mean(1)).mean()
    receptor_mean = np.exp(np.log(receptor_expression).mean(1)).mean()
    interaction_score =  ligand_mean * receptor_mean

    
    # If we specify downstream TF expression, we take the arithmetic mean, as we only need one to be "signficantly" activated
    if tf_expression is not None:
            
        tf_mean = tf_expression.mean()
        interaction_score =  ligand_mean * (receptor_mean * tf_mean) ** 0.5
        
    return interaction_score
