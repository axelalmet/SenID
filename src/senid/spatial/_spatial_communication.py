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
                                sample_key: str = None,
                                min_spots=1,
                                min_spot_pct: float = 0.05,
                                distance_thr: float = 20.0,
                                **kwargs) -> Dict[str, AnnData] | None:
    
    # Calculate the graph-based distance for greater computational efficiency,
    # so we don't have to use Euclidean distance
    if 'spatial_distance' not in adata.obsp:
        _calculate_spatial_neighbors_graph(adata, **kwargs)

    if ligand_receptor_db is None:
        lr_db = ct.pp.ligand_receptor_database(species=model, signaling_type=signaling_type)

    samples = adata.obs[sample_key].unique().tolist() if sample_key is not None else None

    if samples is not None:
        
        output_by_samples = {}
        for sample in samples:
            adata_sample = adata[adata.obs[sample_key] == sample]

            sc.pp.filter_genes(adata_sample, min_cells=min_spots)

            # Filter the cellchat dataframe
            lr_db_filtered = ct.pp.filter_lr_database(lr_db, adata_sample, min_cell_pct=min_spot_pct)

            # Run COMMOT
            ct.tl.spatial_communication(adata_sample,
                                        database_name='cellchat' if ligand_receptor_db is None else db_name,
                                        df_ligrec=lr_db_filtered,
                                        dis_thr=distance_thr,
                                        heteromeric=True,
                                        pathway_sum=True)
        
            output_by_samples[sample] = adata_sample
        
        return output_by_samples
    
    else: # Otherwise there's only one sample

            sc.pp.filter_genes(adata, min_cells=min_spots)

            # Filter the cellchat dataframe
            lr_db_filtered = ct.pp.filter_lr_database(lr_db, adata, min_cell_pct=min_spot_pct)

            # Run COMMOT
            ct.tl.spatial_communication(adata_sample,
                                        database_name='cellchat' if ligand_receptor_db is None else db_name,
                                        df_ligrec=lr_db_filtered,
                                        dis_thr=distance_thr,
                                        heteromeric=True,
                                        pathway_sum=True)