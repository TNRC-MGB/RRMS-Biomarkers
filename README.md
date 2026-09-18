![](docs/cover_image.jpg)

<sub>Cover art by Sarah C. Evans.</sub>

# EBV reactivation priming of the peripheral immune system in multiple sclerosis relapse

Devin A. King, Shrishti Saxena, Danielle Caefer, Kyle C. Downer, Laura E. Saucier,
Ethan Goodman, Jonmichael Aracena, Alena Zhirova, Anthilia Alchanat, Saoirse Nolan,
Hrishikesh Lokhande, Howard L. Weiner, Benjamin E. Gewurz &
Tanuja Chitnis [✉](mailto:tchitnis@bwh.harvard.edu)

*Nature Medicine* (2026) | [Cite the article](https://doi.org/10.1038/s41591-026-04665-3#citeas)

___

**Contact**: [Devin A. King, PhD](mailto:devin.king.neuro@gmail.com) | [ORCID](https://orcid.org/0009-0005-6485-2362)

**scRNA-seq raw sequencing data**: [GSE347634](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE347634)

**CD19+ B cell raw sequencing data**: [GSE344578](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE344578)

**Analysis code**: [*This GitHub repository*](https://github.com/TNRC-MGB/RRMS-Biomarkers)

**Data files**: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22806570.svg)](https://doi.org/10.5281/zenodo.22806570)

Please consider *[opening an issue](https://github.com/TNRC-MGB/RRMS-Biomarkers/issues)*
rather than emailing, to create a public record of how these analyses can be
validated and extended for the benefit of the MS community.

___

## Overview

This repo contains the analysis and figure code for our observational study of
longitudinal blood samples from patients with relapsing-remitting multiple
sclerosis (MS), combining scRNA-seq, CD19+ bulk RNA-seq, flow cytometry and EBV
RT-qPCR into a time-resolved atlas of immune perturbations surrounding relapse.

## Repository layout

```text
R/                           Figure and analysis scripts
bulk/                        CD19+ bulk RNA-seq
scdrs/                       MAGMA and scDRS pipeline
data/                        Small inputs, including gene sets and frozen tables
Supplementary Information/   Supplementary Tables S1 to S9
docs/                        README images
```

The scripts also use additional data and source files deposited in Zenodo
([10.5281/zenodo.22806570](https://doi.org/10.5281/zenodo.22806570)).

## Data availability

| Data | Repository | Accession |
| --- | --- | --- |
| scRNA-seq sequencing data | GEO | [GSE347634](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE347634) |
| CD19+ B cell bulk RNA-seq sequencing data | GEO | [GSE344578](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE344578) |
| Whole-genome sequencing data (WGS) | controlled access | controlled access |
| Genotype-based demultiplexing calls from WGS ([Vireo](https://doi.org/10.1186/s13059-019-1865-2) `donor_ids.tsv`) | Zenodo | [10.5281/zenodo.22806570](https://doi.org/10.5281/zenodo.22806570) |
| CD19+ B cell flow cytometry | Zenodo | [10.5281/zenodo.22806570](https://doi.org/10.5281/zenodo.22806570) |
| qPCR | Zenodo | [10.5281/zenodo.22806570](https://doi.org/10.5281/zenodo.22806570) |
| Transcript quantification/counts (Salmon and Cell Ranger) | Zenodo | [10.5281/zenodo.22806570](https://doi.org/10.5281/zenodo.22806570) |
| Additional data files (scDist fits, NEBULA table, etc) | Zenodo | [10.5281/zenodo.22806570](https://doi.org/10.5281/zenodo.22806570) |

All analyses that depend on whole-genome sequencing can be fully reproduced
from the Vireo donor ID calls provided in Zenodo. All data are fully
de-identified. MS GWAS summary statistics are third-party; see
[`scdrs/README.md`](scdrs/README.md).

## Software environment

R v4.5.2. Shell pipelines in `bulk/` and `scdrs/` were run under WSL Ubuntu
24.04.3 LTS with command-line tools in a conda environment named `rrms`. Every
script that uses the random number generator sets a seed. However, please note
that several steps are sensitive to package versions.

## Running the pipeline

1. **Bulk RNA-seq**, [`bulk/README.md`](bulk/README.md): decoy-aware Salmon
   index on `refdata-gex-GRCh38-2020-A`, 130 CD19+ libraries, TMM-normalized
   `DGEList`.
2. **scRNA-seq preprocessing**: `R/preprocessing_pbmc.R`, then
   `R/preprocessing_bcells.R`. Writes `objects/pbmc_final.RDS` (about 5.5 GB)
   and the B-cell object.
3. **scDist**: `R/scrna_pbmc_scdist.R` writes the five whole-PBMC fits to
   `zenodo/scdist/`; `R/preprocessing_bcells.R` writes the B-cell fits.
4. **scDRS**, [`scdrs/README.md`](scdrs/README.md): MAGMA gene analysis on the
   MS GWAS, then per-cell disease-relevance scores.
5. **Figures**: the scripts in [`R/README.md`](R/README.md), saved to
   `Intermediate/`.

## Figures

| Figure | Script | Panels |
| --- | --- | --- |
| 1. Study design | [`R/figure1.R`](R/figure1.R) | **a** cohort and assay overview; **b** PBMC samples by disease state; **c** age, sex and disease-modifying therapy distributions; **d** temporal sampling design relative to relapse onset; **e** assay coverage matrix |
| 2. Time-resolved single-cell atlas | [`R/figure2.R`](R/figure2.R) | **a** PBMC UMAP by cell state; **b** marker heatmaps supporting annotation; **c** scDist perturbation UMAPs, three contrasts; **d** cell-state scDist distances; **e** patient-aware scDRS on MS GWAS; **f** UMAP of mean scDRS; **g** scDist distance vs pre-relapse scDRS change |
| 3. Atypical B-cell axis | [`R/figure3.R`](R/figure3.R) | **a** UMAP of 14 re-clustered B-cell subsets; **b** RNA marker dot plot; **c** imputed surface-protein overlays; **d** subset-level scDist distances; **e** patient-aware scDRS across subsets; **f** scDist distance vs scDRS change |
| 4. ABC-like and gp350+ expansion | [`R/figure4.R`](R/figure4.R) | **a** flow cytometry workflow schematic, not in R; **b** FlowSOM metacluster UMAP; **c** core and overlay marker heatmaps; **d** within-patient ABC/ABC-like log2 fold change, 23 pairs; **e** gp350+ cells on UMAP; **f** gp350+ frequency per metacluster; **g** gp350+ frequency by disease state |
| 5. Viral and host programs | [`R/figure5.R`](R/figure5.R), [`bulk/figure5e.R`](bulk/figure5e.R) | **a** 40 - dCt per EBV transcript; **b** PCA on EBV transcripts; **c** odds ratios for imminent relapse; **d** over-representation of the pre-relapse ABC scDist signature in EBV-host modules; **e** bulk LMP1 host-response signature over time to relapse |
| 6. Genetic susceptibility and EBV | [`R/figure6.R`](R/figure6.R) | **a** EBNA2 and KSHV LANA NES heatmap across B-cell subsets; **b** scDist distance vs mean EBNA2 NES, colored by scDRS; **c** working model, BioRender |
| Supplementary 1. Sample utilization | [`R/supplementary_figure1.R`](R/supplementary_figure1.R) | **a** sample collection times relative to relapse; **b** UpSet of assay intersections; **c** samples per patient |
| Extended Data 1 to 10 | `R/extended_data_*.R` | see [`R/README.md`](R/README.md) |

## License

Source code is MIT licensed; see [`LICENSE`](LICENSE). Data and third-party
materials keep their own terms.

## References

1. Zhao, B., Zou, J., Wang, H., Johannsen, E., Peng, C., Quackenbush, J., …
   Kieff, E. Epstein-Barr virus exploits intrinsic B-lymphocyte transcription
   programs to achieve immortal cell growth. Proc. Natl. Acad. Sci. U.S.A. 108,
   14902–14907 (2011). https://doi.org/10.1073/pnas.1108892108
2. Arvey, A., Tempera, I., Tsai, K., Chen, H.-S., Tikhmyanova, N., Klichinsky,
   M., … Lieberman, P. M. An Atlas of the Epstein-Barr Virus Transcriptome and
   Epigenome Reveals Host-Virus Regulatory Interactions. Cell Host & Microbe 12,
   233–245 (2012). https://doi.org/10.1016/j.chom.2012.06.008
3. Mercier, A., Arias, C., Madrid, A. S., Holdorf, M. M. & Ganem, D.
   Site-Specific Association with Host and Viral Chromatin by Kaposi's
   Sarcoma-Associated Herpesvirus LANA and Its Reversal during Lytic
   Reactivation. J Virol 88, 6762–6777 (2014).
   https://doi.org/10.1128/JVI.00268-14
4. De Leeuw, C. A., Mooij, J. M., Heskes, T. & Posthuma, D. MAGMA: Generalized
   Gene-Set Analysis of GWAS Data. PLoS Comput Biol 11, e1004219 (2015).
   https://doi.org/10.1371/journal.pcbi.1004219
5. Van Gassen, S., Callebaut, B., Van Helden, M. J., Lambrecht, B. N.,
   Demeester, P., Dhaene, T. & Saeys, Y. FlowSOM: Using self‐organizing maps for
   visualization and interpretation of cytometry data. Cytometry Pt A 87,
   636–645 (2015). https://doi.org/10.1002/cyto.a.22625
6. Korotkevich, G., Sukhov, V., Budin, N., Shpak, B., Artyomov, M. N. &
   Sergushichev, A. Fast gene set enrichment analysis. Preprint at
   https://doi.org/10.1101/060012 (2016).
7. Patro, R., Duggal, G., Love, M. I., Irizarry, R. A. & Kingsford, C. Salmon
   provides fast and bias-aware quantification of transcript expression. Nat
   Methods 14, 417–419 (2017). https://doi.org/10.1038/nmeth.4197
8. International Multiple Sclerosis Genetics Consortium, Patsopoulos, N. A.,
   Baranzini, S. E., Santaniello, A., Shoostari, P., Cotsapas, C., … De Jager,
   P. L. Multiple sclerosis genomic map implicates peripheral immune cells and
   microglia in susceptibility. Science 365, eaav7188 (2019).
   https://doi.org/10.1126/science.aav7188
9. Korsunsky, I., Millard, N., Fan, J., Slowikowski, K., Zhang, F., Wei, K., …
   Raychaudhuri, S. Fast, sensitive and accurate integration of single-cell data
   with Harmony. Nat Methods 16, 1289–1296 (2019).
   https://doi.org/10.1038/s41592-019-0619-0
10. He, L., Davila-Velderrain, J., Sumida, T. S., Hafler, D. A., Kellis, M. &
    Kulminski, A. M. NEBULA is a fast negative binomial mixed model for
    differential or co-expression analysis of large-scale multi-subject
    single-cell data. Commun Biol 4, 629 (2021).
    https://doi.org/10.1038/s42003-021-02146-6
11. Ashhurst, T. M., Marsh‐Wakefield, F., Putri, G. H., Spiteri, A. G., Shinko,
    D., Read, M. N., … King, N. J. C. Integration, exploration, and analysis of
    high‐dimensional single‐cell cytometry data using Spectre. Cytometry Pt A
    101, 237–253 (2022). https://doi.org/10.1002/cyto.a.24350
12. Zhang, M. J., Hou, K., Dey, K. K., Sakaue, S., Jagadeesh, K. A., Weinand,
    K., … Price, A. L. Polygenic enrichment distinguishes disease associations
    of individual cells in single-cell RNA-seq data. Nat Genet 54, 1572–1580
    (2022). https://doi.org/10.1038/s41588-022-01167-z
13. Viel, K. C. M. F., Parameswaran, S., Donmez, O. A., Forney, C. R., Hass, M.
    R., Yin, C., … Weirauch, M. T. Shared and distinct interactions of type 1
    and type 2 Epstein-Barr Nuclear Antigen 2 with the human genome. BMC
    Genomics 25, 273 (2024). https://doi.org/10.1186/s12864-024-10183-8
14. Nicol, P. B., Paulson, D., Qian, G., Liu, X. S., Irizarry, R. & Sahu, A. D.
    Robust identification of perturbed cell types in single-cell RNA-seq data.
    Nat Commun 15, 7610 (2024). https://doi.org/10.1038/s41467-024-51649-3
15. Zalewski, D. & Bogucka-Kocka, A. RQdeltaCT: an open-source R package for
    relative quantification of gene expression using delta Ct methods. Sci Rep
    15, 29762 (2025). https://doi.org/10.1038/s41598-025-11822-0

EBNA2 gene sets derive from ChIP-seq series
[GSE246060](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246060) and
[GSE29498](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE29498), and
the KSHV LANA gene set from
[GSE56144](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE56144). The
processed sets are in `data/genesets/`.
