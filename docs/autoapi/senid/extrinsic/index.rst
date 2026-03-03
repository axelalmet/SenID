senid.extrinsic
===============

.. py:module:: senid.extrinsic


Functions
---------

.. autoapisummary::

   senid.extrinsic.senchat
   senid.extrinsic.calculate_interaction_scores_for_lr_pair
   senid.extrinsic.constrain_senchat_interactions
   senid.extrinsic.subset_for_communication
   senid.extrinsic.run_cnmf
   senid.extrinsic.infer_consensus_modules
   senid.extrinsic.calculate_ligand_enrichment


Package Contents
----------------

.. py:function:: senchat(adata, sender_label, receiver_label = None, sender_groups = None, receiver_groups = None, output_key = 'SenChat_output', model = 'human', use_highly_variable = False, highly_variable_key = None, stringency = 'neither', min_signal_proportion = 0.1, filter_small_groups = False, min_cell_proportion = 0.001, test_permutation = False, n_perms = 100, return_df = False, n_jobs = 1, score_method = 'v1', perm_method = 'v1')

   Calculate all interaction scores with respect to a specified lbel.

   :param adata: Annotated data matrix containing expression data.
   :type adata: AnnData
   :param group_label: Key in `adata.obs` containing the group label.
   :type group_label: str
   :param model: The model organism. Options are 'mouse' or 'human'.
   :type model: str
   :param use_highly_variable: Whether to use highly variable genes. If True, then you need to have previously determined which genes are highly variable, i.e.,
                               the information shoudl be contained in adata.var[highly_variable_key].
   :type use_highly_variable: bool
   :param highly_variable_key: Key in adata.var containing the highly variable genes. If None, then the highly variable genes are not used.
   :type highly_variable_key: Optional[str]
   :param min_proportion: Minimum proportion of cells expressing the ligand or receptor to be considered as a potential sender or receiver.
   :type min_proportion: float
   :param test_permutation: Whether to calculate an associated p-value for the interaction scores using permutation testing.
   :type test_permutation: bool
   :param n_perms: Number of permutations to use for permutation testing
   :type n_perms: int

   :returns: DataFrame containing the interaction scores across all ligand-receptor interaction pairs and potentially interacting groups
   :rtype: pd.DataFrame


.. py:function:: calculate_interaction_scores_for_lr_pair(adata, ligand, receptor, sender_label, receiver_label = None, sender_groups = None, receiver_groups = None, downstream_tf = None, pathway = None, pathway_type = None, is_neurotransmitter = None, min_signal_proportion = 0.1, filter_small_groups = False, min_cell_proportion = 0.001, test_permutation = False, n_perms = 100, score_method = 'v1', perm_method = 'v1')

   Calculate all interaction scores for a specified ligand-receptor pair

   :param adata: Annotated data matrix containing expression data.
   :type adata: AnnData
   :param sender_label: Key in `adata.obs` containing the group label.
   :type sender_label: str
   :param group1: Name of the first group.
   :type group1: str
   :param group2: Name of the second group.
   :type group2: str
   :param ligand_receptor_pairs: DataFrame containing ligand-receptor pairs.
   :type ligand_receptor_pairs: pd.DataFrame
   :param interaction_score_method: Method to calculate the interaction score. Options are "geometric_mean" or "arithmetic_mean".
   :type interaction_score_method: str
   :param tf_expression: Key in `adata.obs` containing the expression of a downstream transcription factor. If None, the interaction score is calculated without considering the downstream transcription factor.
   :type tf_expression: Optional[Union[str, np.ndarray]]

   :returns: DataFrame containing the interaction scores for a specified ligand-receptor pair
             between all potentially interacting groups.
   :rtype: pd.DataFrame


.. py:function:: constrain_senchat_interactions(adata, senchat_output_key, implausible_interactions = None, pathway_types = None)

.. py:function:: subset_for_communication(adata, senchat_output_key, pval_threshold = None, subset_sasp = True, layer = 'counts', model = 'human')

.. py:function:: run_cnmf(output_dir, adata_fname, output_name, n_modules, combine_only = False, worker_i = 0, total_workers = 1, **kwargs)

.. py:function:: infer_consensus_modules(adata_cnmf, output_dir, output_name, optimal_k, density_threshold = 0.1)

.. py:function:: calculate_ligand_enrichment(adata, spectra_tpm_key = 'cNMF_spectra_tpm')

