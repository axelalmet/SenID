from typing import Literal, Dict
from anndata import AnnData
import commot as ct
import scanpy as sc
import pandas as pd
import squidpy as sq

def _calculate_spatial_neighbors_graph(adata: AnnData,
                                       **kwargs) -> None:
     
     # Calculate the spatial graph
     sq.gr.spatial_neighbors(adata, **kwargs)

def infer_spatial_communication(adata: AnnData,
                                ligand_receptor_db: pd.DataFrame = None,
                                db_name: str = None,
                                model: Literal['human', 'mouse'] = 'human',
                                signaling_type: Literal['Secreted Signaling', 'Cell-Cell Contact', 'ECM-Receptor'] = 'Secreted Signaling',
                                filter_spots: bool = None,
                                min_spots: int = 1,
                                min_spot_pct: float = 0.05,
                                distance_thr: float = 20.0,
                                **kwargs) -> None:
    
    # Calculate the graph-based distance for greater computational efficiency,
    # so we don't have to use Euclidean distance
    if 'spatial_distance' not in adata.obsp:
        _calculate_spatial_neighbors_graph(adata, **kwargs)

    if filter_spots: 
        sc.pp.filter_genes(adata, min_cells=min_spots)

    if ligand_receptor_db is None:
        lr_db = ct.pp.ligand_receptor_database(species=model, signaling_type=signaling_type)

        # Filter the cellchat dataframe
        lr_db_filtered = ct.pp.filter_lr_database(lr_db, adata, min_cell_pct=min_spot_pct)

        # Run COMMOT
        ct.tl.spatial_communication(adata,
                                    database_name='cellchat' if ligand_receptor_db is None else db_name,
                                    df_ligrec=lr_db_filtered,
                                    dis_thr=distance_thr,
                                    heteromeric=True,
                                    pathway_sum=True)    

    else: # Otherwise there's only one sample

            # Filter the cellchat dataframe
            lr_db_filtered = ct.pp.filter_lr_database(lr_db, adata, min_cell_pct=min_spot_pct)

            # Run COMMOT
            ct.tl.spatial_communication(adata,
                                        database_name='cellchat' if ligand_receptor_db is None else db_name,
                                        df_ligrec=lr_db_filtered,
                                        dis_thr=distance_thr,
                                        heteromeric=True,
                                        pathway_sum=True)