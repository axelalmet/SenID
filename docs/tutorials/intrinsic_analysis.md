# Intrinsic analysis with SenID

This tutorial demonstrates how to identify intrinsic senescence-associated transcriptional dynamics using Monod inference within SenID.

We use scRNA-seq of damaged lung tissue from 18-month-old mice, originally generated in [Reyes et al. (2022)](https://doi.org/10.1126/science.abf3326). In this study, senescent cells were labelled using GFP markers of p16^INK4a^, a common senescence marker. The authors then performed FACS prior to sequencing to generate GFP+ and GFP- scRNA-seq datasetes.

We re-aligned the data using [kb-python](https://kallisto.readthedocs.io/en/latest/) to generate unspliced and spliced counts, and re-annotated the data based on spliced expression using a procedure similar to the scANVI [tutorial](https://docs.scvi-tools.org/en/stable/tutorials/notebooks/scrna/tabula_muris.html) that uses the Tabula Muris Senis dataset as a reference.

## Load the data

```python
from pathlib import Path
import scanpy as sc
import senid as sid

adata = sc.read_h5ad("reyes22_old_unspliced.h5ad")
```

We will analyze intrinsic dynamic changes between GFP-positive senescent cells and GFP-negative per cell type using [Monod](https://github.com/pachterlab/monod). We therefore need to stratify the data by cell type and GFP status.

## Create output directories

```python
output_dir = Path("intrinsic_outputs")
output_dir.mkdir(exist_ok=True)
```

## Generate loom objects

SenID partitions cells by cell type and condition before preparing loom files for Monod inference. Here, `partition_key` is used to split by cell type, and `comparison_key` is the key to split cell-type-specific populations into GFP+ and GFP- subsets. This step generates a list of all files that are to be used for
Monod inference.


```python
cell_groups = sid.it.generate_loom_objects(
    adata,
    dataset_key="reyess2_old",
    partition_key="prediction_scanvi",
    comparison_key="GFP",
    output_directory=output_dir,
    min_cells=100,
    return_groups=True,
)
```


## Prepare Monod inference files

```python
sid.it.prepare_files_for_inference(
    dataset_key="reyes22_old",
    output_directory=output_dir,
    cell_groups=cell_groups,
)
```

This step prepares input matrices and metadata required for Monod.

---

## Run Monod inference

```python
sid.it.run_monod_inference(
    dataset_key="reyes22_old",
    output_directory=output_dir,
    cell_groups=cell_groups,
)
```

This step estimates transcriptional dynamics for each cell group.

Depending on dataset size, this may require substantial compute time.

---

## Process Monod fits

```python
sid.it.process_monod_fits(
    dataset_key="reyes22_old",
    output_directory=output_dir,
    cell_groups=cell_groups,
)
```

This step extracts interpretable transcriptional dynamics statistics from the raw Monod outputs.

---

## Identify differential transcriptional dynamics

```python
de_results = sid.it.get_monod_de_genes(
    dataset_key="reyes22_old",
    output_directory=output_dir,
    comparison_1="GFPPos",
    comparison_2="GFPNeg",
)
```

This function identifies genes exhibiting differential transcriptional dynamics between conditions.

---

## Save outputs

```python
de_results.to_csv(
    output_dir / "intrinsic_differential_dynamics.csv",
    index=False,
)
```