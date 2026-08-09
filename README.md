# Reanalysis of a glycolysis-weighted expression contrast in human ICH

This repository reproduces the evidence-aligned secondary analysis of longitudinal sorted CD14+ monocyte/macrophage transcriptomes from **GSE163256**. It corrects a prior data-provenance error: GSE163256 is sorted-cell transcriptomics, not single-cell RNA sequencing. **GSE166638** is a one-patient scRNA-seq case study and is not used as an independent validation cohort.

## Primary analysis

- Input: `GSE163256_monos_log_fpkm_techreps_collapsed.csv.gz`
- Biological unit: patient
- Locked ICH subset: 136 blood/hematoma observations from 20 patients, days 1-6
- Duplicate handling: three residual duplicate sample labels are collapsed by their mean
- Predictor: mean Hallmark glycolysis expression minus mean Hallmark oxidative-phosphorylation expression
- Outcomes: mean Hallmark inflammatory-response and hypoxia expression
- Primary model: `outcome ~ contrast + compartment + day + (1 | patient)`
- Sensitivity checks: patient-clustered robust covariance, cross-signature overlap pruning, glycolysis/OXPHOS component decomposition, leave-one-gene-out refitting, and leave-one-patient-out prediction

The expression contrast is an exploratory transcriptomic score. It does not measure metabolic flux.

## Reproduce

Requirements: Conda/Mamba and Snakemake 9 or later; R is required only for the summary figure.

```bash
conda env create -f environment.yml
conda activate ich-expression-contrast
snakemake -n --snakefile workflow/Snakefile --configfile config/config.yaml
snakemake --cores 1 --snakefile workflow/Snakefile --configfile config/config.yaml
```

The workflow downloads the public GEO matrix and MSigDB Hallmark GMT release specified in `config/config.yaml`. If the MSigDB download endpoint requires acceptance of updated terms, download the named GMT manually and place it at `data/raw/h.all.v2025.1.Hs.symbols.gmt`.

## Clinical-variable availability

The GSE163256 cohort came from the MISTIE III surgical arm. Catheter placement, hematoma drainage, and repeated rtPA administration were design features, and the source study explicitly compared samples collected before and after rtPA exposure. Public expression-matrix labels do not provide observation-linked antibiotic, corticosteroid, infection, or detailed medication histories. These variables therefore cannot be added reliably to the present model. See `provenance/clinical_variable_audit.md`.

## Repository contents

- `workflow/Snakefile`: complete rule graph
- `scripts/`: download, scoring, model-fitting, and R figure scripts
- `config/config.yaml`: public URLs and versioned parameters
- `results/`: frozen numerical outputs
- `figures/`: publication summary figure
- `provenance/`: source and clinical-variable audit

`provenance/reference_verification.csv` records the 40 manuscript references, DOI or official URL, verification link, status, and live-check date (2026-08-09).

## Data sources

- GSE163256: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163256
- GSE166638: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166638
- Source report for GSE163256: https://doi.org/10.1126/sciimmunol.abd6279
- Source report for GSE166638: https://doi.org/10.1172/jci.insight.145857

## License and attribution

Analysis code is released under the MIT License. Source expression data remain governed by GEO/source-study terms. MSigDB content is not redistributed here; users should comply with the current MSigDB license and cite the Hallmark collection.
