senid.intrinsic
===============

.. py:module:: senid.intrinsic


Functions
---------

.. autoapisummary::

   senid.intrinsic.generate_loom_objects
   senid.intrinsic.prepare_files_for_inference
   senid.intrinsic.run_monod_inference
   senid.intrinsic.process_monod_fits
   senid.intrinsic.get_monod_de_genes


Package Contents
----------------

.. py:function:: generate_loom_objects(adata, dataset_key, partition_key, comparison_key = None, str_replace = '-', output_directory = '.', min_cells = 100, return_groups = False, **kwargs)

   Generate loom files for each phenotype in the AnnData object.

   Parameters:
   - adata: AnnData object containing the data.
   - pheno_col: Column name in adata.obs that contains phenotype information.
   - data_directory: Directory to save the loom files.


.. py:function:: prepare_files_for_inference(loom_directory, dataset_names, transcriptome_filepath)

.. py:function:: run_monod_inference(loom_directory, dataset_names, transcriptome_filepath, dir_string, dataset_strings, phys_lb, phys_ub, samp_lb, samp_ub, gridsize, max_iterations = 15, n_restarts = 1, n_jobs = 1)

.. py:function:: process_monod_fits(sr_arr, sd_arr, dir_string, plot_results = True, n_jobs = 1)

.. py:function:: get_monod_de_genes(sr1, sr2, sd1, sd2, dname1, dname2, gene_names, gf_rej = False, param_lfc = 2.0, mean_lfc = 1.0, pval_thr = 0.05, outlier_de = True, single_nuc = False, correct_off = False)

