# CD19+ bulk RNA-seq quantification

Quantifies 130 CD19+ bulk RNA-seq libraries against a decoy-aware Salmon index
and imports them to a TMM-normalised `DGEList`.

The shell stages were run under WSL Ubuntu 24.04.3 LTS with the conda 
environment `rrms`.

## Code

| Script | Purpose | Inputs | Key Outputs |
| --- | --- | --- | --- |
| `01_build_index.sh` | builds a decoy-aware Salmon index | `refdata-gex-GRCh38-2020-A` genome and GTF | `bulk/reference/salmon_index/idx/`, `versions.txt` |
| `02_quantify.sh` | quantifies each sample with Salmon | `bulk/fastq/`, the index | `zenodo/bulk/counts/<Sample_ID>/quant.sf`, `manifest.tsv`, `02_quantify.log` |
| `import_counts.R` | tximport, filtering, TMM normalization | `zenodo/bulk/counts/`, `genes.gtf` | `bulk/results/dge.rds`, `txi.rds`, `gene_counts.tsv` |
| `figure5e.R` | renders Figure 5e | `bulk/results/dge.rds`, S7, `pebv500.RDS` | `figure5e_bigpoints.pdf`, `figure5e_scores.tsv`, `figure5e_qc_metrics.tsv` |

```bash
bash bulk/01_build_index.sh
JOBS=3 THREADS=6 bash bulk/02_quantify.sh
Rscript bulk/import_counts.R
Rscript bulk/figure5e.R
```

## Inputs

FASTQs are in GEO. Transcript quantification results are in Zenodo.

- `bulk/fastq/` - 260 files, paired-end, named
  `<Sample_ID>_S*_L*_R{1,2}_001.fastq.gz`. The `Sample_ID` prefix corresponds to
  Supplementary Tables S3 and S7.
- `bulk/reference/refdata-gex-GRCh38-2020-A/`. The 10x GRCh38 reference.
- `Supplementary Information/S7 B Cell Bulk RNA-seq Libraries.xlsx`, listing the
  130 libraries.

## Outputs

| Path | Contents | In this repo? |
| --- | --- | --- |
| `bulk/reference/salmon_index/idx/` | Salmon index | no |
| `bulk/reference/salmon_index/versions.txt` | Salmon and gffread versions, reference, k-mer, build date | no |
| `zenodo/bulk/counts/<Sample_ID>/quant.sf` | per-sample quantification | Zenodo |
| `bulk/results/dge.rds` | filtered, TMM-normalized `DGEList` | no |
| `bulk/results/gene_counts.tsv` | gene-level count matrix | no |
| `bulk/results/figure5e_*.pdf`, `*.tsv` | Figure 5e and its source values | no |

