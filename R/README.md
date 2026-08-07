# Figures

| Figure | File | Panels |
| --- | --- | --- |
| Figure 1 | `figure1.R` | a-e, cohort |
| Figure 2 | `figure2.R` | a-g, PBMC atlas, scDist, scDRS |
| Figure 3 | `figure3.R` | a-f, B-cell compartment |
| Figure 4 | `figure4.R` | b-g, FlowSOM via `facs_flowsom_pipeline.R` |
| Figure 5 | `figure5.R` | a-d, qPCR and scDist ORA; see `bulk/` for e |
| Figure 6 | `figure6.R` | a-b, EBV/EBNA2 |
| Supplementary 1 | `supplementary_figure1.R` | a-c, cohort supplement |
| Extended Data 1 | `extended_data_1.R` | scRNA QC pipeline |
| Extended Data 2 | `extended_data_2.R` | Azimuth annotation |
| Extended Data 3 | `extended_data_3.R` | scRNA cohort composition |
| Extended Data 4 | `extended_data_4.R` | scDist null / permutation |
| Extended Data 5 | `extended_data_5.R` | scDist relapse and pre-relapse dotplots |
| Extended Data 6 | `extended_data_6.R` | B-cell subclustering / QC |
| Extended Data 7 | `extended_data_7.R` | B-cell DEG volcano grid |
| Extended Data 8 | `extended_data_8.R` | flow supplement, b-e |
| Extended Data 9 | `extended_data_9.R` | per-subset EBV GSEA dot plots |
| Extended Data 10 | `extended_data_10.R` | NES by EBV lifecycle stage |

`facs_flowsom_pipeline.R` is the FlowSOM pipeline for Figure 4 and
Extended Data 8

## Run order

1. **Pipelines** build the shared objects the figure scripts load e.g.: These 
objects can also be downloaded from Zenodo or requested directly from us.
   - `preprocessing_pbmc.R` writes `objects/pbmc_final.RDS`
   - `scrna_pbmc_scdist.R` writes the five whole-PBMC scDist fits to
     `zenodo/scdist/`
   - `preprocessing_bcells.R` writes `objects/B_clean_scdrs.RDS` and the
     B-cell scDist fits
2. **Figure scripts** generate the plots from the objects above.

## Reproducibility

Every script that uses the RNG sets a seed. Three things are seeded but still 
sensitive to the package version, so they are worth checking after any package 
upgrade:

- `facs_flowsom_pipeline.R` trains a SOM iteratively in floating point, and
  `metaClustering_consensus` can flip a near-tie between metaclusters.
- `fgseaMultilevel`, used by Figure 6 and Extended Data 9 and 10, estimates
  p-values by sampling.
- see relevant notes for scDblFinder in scRNA pipeline

<!--
```r
source("R/00_environment.R")
check_environment(); pin_environment(); record_versions()
```
-->