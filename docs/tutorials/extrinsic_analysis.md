# Extrinsic communication analysis with SenID
This tutorial demonstrates how to identify extrinsic senescence-associated patterns using SenID.

We use scRNA-seq of damaged lung tissue from 18-month-old mice, originally generated in [Reyes et al. (2022)](https://doi.org/10.1126/science.abf3326). In this study, senescent cells were labelled using GFP markers of p16^INK4a^, a common senescence marker. The authors then performed FACS prior to sequencing to generate GFP+ and GFP- scRNA-seq datasetes.

We re-aligned the data using [kb-python](https://kallisto.readthedocs.io/en/latest/) to generate unspliced and spliced counts, and re-annotated the data based on spliced expression using a procedure similar to the scANVI [tutorial](https://docs.scvi-tools.org/en/stable/tutorials/notebooks/scrna/tabula_muris.html) that uses the Tabula Muris Senis dataset as a reference. To infer cell-cell communication, we normalised spliced counts by estimated library sizes using [scran](https://bioconductor.org/packages/release/bioc/html/scran.html) and then log-transformed the normalised counts as a variance-stabilising transform.

## Load the data

```python
from pathlib import Path
import scanpy as sc
import senid as sid

adata = sc.read_h5ad("reyes22_old_unspliced.h5ad")
```

We will analyze extrinsic communication patterns. 

## Create output directories

```python
output_dir = Path("extrinsic_outputs")
output_dir.mkdir(exist_ok=True)
```

## Select highly deviant genes for SenChat

Before calculating communication scores, we select highly deviant genes from the spliced expression layer.

```python
sid.pp.binomial_deviance_selection(
    adata,
    layer="spliced",
    deviance_key="binomial_deviance",
    highly_variable_key="highly_deviant",
    n_top_genes=2000,
)
```

## Infer cell-cell communication using SenChat

SenChat calculates ligand-receptor interaction scores between sender and receiver cell populations. Here, both sender and receiver populations are defined using the same cell type annotation column, `predictions_scanvi`.

```python

sid.tl.calculate_interaction_scores(
    adata,
    sender_label="predictions_scanvi",
    receiver_label="predictions_scanvi",
    model="mouse",
    use_highly_variable=True,
    highly_variable_key="highly_deviant",
    stringency="neither",
    min_proportion=0.2,
    test_permutation=True,
    n_perms=100,
    return_df=False,
    n_jobs=4,
)
```

The output will be saved into a `.uns` key, `senchat_output_key = "SenChat_output_predictions_scanvi_predictions_scanvi"`, that we will need for communication analysis.

```python
senchat_output_key = "SenChat_output_predictions_scanvi_predictions_scanvi"
```

SenID can remove biologically implausible interactions and restrict the analysis to selected communication classes. The former are used to remove cell types that we know aren't close enough, e.g., in tissues like retina, which are highly structured, certain cell types will definitely not be in contact with others.

```python
sid.et.constrain_senchat_interactions(
    adata,
    senchat_output_key=senchat_output_key,
    implausible_interactions=None,
    pathway_types=[
        "Secreted Signaling",
        "Cell-Cell Contact",
    ],
)
```

Here, we retain secreted signalling and cell-cell contact interactions.

```python
constrained_senchat_output_key = f"{senchat_output_key}_constrained"
```

## Subset for senescence-associated communication genes

Next, we subset the AnnData object to genes involved in significant communication interactions. Setting `subset_sasp=True` restricts the analysis to senescence-associated secretory phenotype (SASP)-related communication genes.

```python
adata_ccc = sid.et.subset_for_communication(
    adata,
    senchat_output_key=constrained_senchat_output_key,
    pval_threshold=0.05,
    subset_sasp=True,
    layer="spliced",
    model="mouse",
)
```

Save this subset for future use.
```python
adata_ccc.write(
    output_dir / "reyes22_old_ccc_sasp_SenID.h5ad",
    compression="gzip",
)
```
## Learn communication modules using cNMF

SenID uses [cNMF](https://github.com/dylkot/cNMF) to decompose communication-associated gene expression into communication programs.

```python
sid.et.run_cnmf(
    output_dir=output_dir,
    adata_fname="reyes22_old_ccc_sasp_SenID.h5ad",
    output_name="reyes22_old_ccc_sasp_SenID",
    n_modules=range(2, 11),
    combine_only=False,
    worker_i=0,
    total_workers=1,
)
```
## Infer consensus communication programs
There will be a file, `{output_name}.k_selection.png`, that indicates what the optimal factorisation choice is, based on stability and error. In this case, the optimal choice seems to be `k=6`. Based on this, we then generate the "consensus" factorisation for `k = 6`, and set the density threshold low enough to remove any clear outliers. You can play with this parameter based on the results shown in `{output_name}.clustering.k_{optimal_k}.dt_{density_threshold.replace(".", "_")}.png`

```python
sid.et.infer_consensus_modules(
    adata_cnmf=adata_ccc,
    output_dir=output_dir,
    output_name="reyes22_old_ccc_sasp_SenID",
    optimal_k=6,
    density_threshold=0.075,
)
```

Save the output.
```python
adata_ccc.write(
    output_dir / "reyes22_old_ccc_sasp_SenID.h5ad",
    compression="gzip",
)
```