senid.preprocessing
===================

.. py:module:: senid.preprocessing


Functions
---------

.. autoapisummary::

   senid.preprocessing.binomial_deviance_selection
   senid.preprocessing.calculate_binomial_deviance_batch


Package Contents
----------------

.. py:function:: binomial_deviance_selection(adata, layer = None, deviance_key = 'binomial_deviance', highly_variable_key = 'highly_deviant', n_top_genes = 1000, batch_key = None, sort_genes = True)

   Python implementation of the brilliantly effective feature selection method, developed by Will Townes
   (see Townes et al. 2019: doi.org/10.1186/s13059-019-1861-6). The idea is that we use a binomial
   deviance to quantify the variability of a gene, based on a multinomial model of UMI counts. We only
   calculate the binomial deviance, whereas the scry package developed by Townes has the option to calculate
   the Poisson deviance.

   :param adata: Annotated data matrix.
   :type adata: AnnData
   :param highly_vairable_key: The key in adata.var to store the highly variable genes.
   :type highly_vairable_key: str, optional (default: 'highly_deviant')
   :param layer: The layer of the AnnData object to use. If None, this method won't work.
   :type layer: str, optional (default: None)
   :param n_top_genes: Number of top genes to select when assigning genes as 'highly variable'
   :type n_top_genes: int, optional (default: 1000)
   :param batch_key: The batch label in adata.obs. If used, we calculate the binomial deviance per batch and
                     then define the binomial deviance per gene as the sum of the per-batch deviances.
   :type batch_key: str, optional (default: None)
   :param sort_genes: If True, sort genes by binomial deviance.
   :type sort_genes: bool, optional (default: True)

   :returns: **adata** -- Annotated data matrix with the binomial deviance per gene stored in adata.var['binomial_deviance'].
             The top_n_genes highly variable genes are stored in adata.var[highly_variable_key].
   :rtype: AnnData


.. py:function:: calculate_binomial_deviance_batch(counts, size_factors, batch_keys = None)

