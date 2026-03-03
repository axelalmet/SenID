senid.spatial
=============

.. py:module:: senid.spatial


Functions
---------

.. autoapisummary::

   senid.spatial.infer_spatial_communication
   senid.spatial.subset_for_spatial_sasp_communication


Package Contents
----------------

.. py:function:: infer_spatial_communication(adata, ligand_receptor_db = None, db_name = None, model = 'human', signaling_type = 'Secreted Signaling', filter_spots = None, min_spots = 1, min_spot_pct = 0.05, distance_thr = 20.0, **kwargs)

.. py:function:: subset_for_spatial_sasp_communication(adata, layer = 'counts', db_name = None, model = 'human')

   Subset the AnnData object for spatial SASP communication analysis.


