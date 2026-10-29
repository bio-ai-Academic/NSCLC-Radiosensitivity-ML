# NSCLC radiosensitivity: DNA repair and oxidative stress

Analysis code associated with **A dual-axis DNA repair and oxidative stress network governs radiosensitivity in non-small cell lung cancer**.

## Scope

The supplied scripts implement transcriptomic data preparation, nested cross-validation for LASSO regression, multilayer perceptron training with DeepSHAP, candidate-gene integration, WGCNA, differential expression, functional enrichment, and sensitivity analyses.

This repository contains code only. Input datasets, trained models, manuscript files, and reported results are not included. The workflow has not been executed on the study data as part of this repository preparation; exact reproduction of manuscript values has not been verified. Data acquisition and ComBat batch correction are not implemented in these scripts: the expression input must already be normalized and appropriately processed.

## Dependencies

Install Python packages using:

```bash
python -m pip install -r requirements.txt
```

R packages: `WGCNA`, `limma`, `clusterProfiler`, `org.Hs.eg.db`, and `enrichplot` (and their dependencies). Bioconductor packages can be installed using `BiocManager::install()`.

Dependency versions are not pinned in the supplied package. TensorFlow/SHAP compatibility should be checked in the target environment.

## Required inputs

Create an `analysis_inputs/` folder next to the scripts:

| File | Expected format |
| --- | --- |
| `expression_log2.csv` | Genes in rows, samples in columns; first column contains gene identifiers. Finite normalized log2 expression values. Use human gene symbols for downstream enrichment. |
| `radiation_response.csv` | Columns `sample_id`, `SF2`, `AUC_recomputed`; optional `AUC_published`, `alpha`, `beta`. Repeated radiation measurements are aggregated by median. |
| `sample_metadata.csv` | Columns `sample_id`, `histology`; IDs must match expression and response data. |

At least 15 matched samples are required by the preparation script. The manuscript describes a RadioGx Cleveland dataset analysis of 54 NSCLC cell lines; the data are not bundled here.

Optional external association inputs are `external_expression_log2.csv` (genes by samples) and `external_response.csv` (`sample_id`, `SF2`) in the same folder.

## Run

From PowerShell, with Python and Rscript available on PATH:

```powershell
.\10.run_workflow.ps1
# Include permutation, dimension/endpoint, external-association and PCA/KernelSHAP sensitivity analyses:
.\10.run_workflow.ps1 -Extended
```

Alternatively, run scripts 1 through 7 in numerical order with `Rscript` for R files and `python` for Python files, then optionally run scripts 8 and 9. Keep `workflow.py` alongside the Python scripts.

## Script guide

| Script | Purpose |
| --- | --- |
| `1.prepare_transcriptomic_data.R` | Match samples, filter expression, summarize measured response, define SF2 groups, and generate QC outputs. |
| `2.lasso_nested_cross_validation.py` | Nested LASSO tuning, out-of-fold prediction, stability resampling; optional `--loocv`. |
| `3.mlp_training_and_deepshap.py` | Five-fold MLP evaluation, inner holdout hyperparameter search, and held-out DeepSHAP attribution. |
| `4.candidate_gene_integration.py` | Fuse LASSO stability and SHAP importance, rank candidates, assess overlap and correlated subsets. |
| `5.wgcna_network_analysis.R` | Co-expression modules and associations with radiation-response traits. |
| `6.differential_expression.R` | limma analyses using SF2-based groups. |
| `7.functional_enrichment.R` | Human gene-set enrichment with GO/KEGG outputs. |
| `8.robustness_and_external_associations.py` | LASSO dimension/endpoint checks, label permutation, and optional external gene associations. |
| `9.pca_kernelshap_sensitivity.py` | PCA reconstruction and gene-level KernelSHAP sensitivity analysis. |
| `10.run_workflow.ps1` | Sequential workflow launcher. |
| `workflow.py` | Shared preprocessing, model, evaluation, and reporting functions. |

## Outputs and interpretation

Results are written under `analysis_outputs/`, including prepared matrices, fold metrics, out-of-fold predictions, model files, gene rankings, enrichment tables, and diagnostic PDFs. Random seed defaults to 123. Results depend on the supplied data and software environment; repository documentation does not certify the manuscript's numerical claims. Cell-line associations require independent validation before clinical application.
