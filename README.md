# 5prime-scrna Pipeline

**Multi-modal analysis of 5' single-cell RNA-seq data from B cells**

A Nextflow pipeline for integrated analysis of 10x Chromium 5' scRNA-seq data, combining gene expression, B cell receptor repertoire, transcription start site, and enhancer analyses. Designed for CVID vs healthy control B cell research.

---

## Overview

This pipeline takes raw FASTQ files from 10x Chromium 5' scRNA-seq experiments (GEX + BCR libraries) and runs the following analyses in parallel:

```
Samplesheet (sample, lane, condition, gex_r1, gex_r2, bcr_r1, bcr_r2)
│
├── FastQC (GEX + BCR) ──────────────────────────────────► MultiQC report
│
├── GEX library
│   ├── Cell Ranger count ──► BAM ──┬── SCAFE ───────────► tCRE matrix (TSS)
│   │                               ├── CamoTSS ──────────► alternative TSS
│   │                               └── GEX matrix
│   │                                       └── Seurat ───► clustering + UMAP
│   └── ReapTEC (STARsolo) ──────────────────────────────► btcEnh → DESeq2
│
└── BCR library
    ├── Cell Ranger VDJ ─────────────────────────────────► clonotypes
    └── MiXCR ───────────────────────────────────────────► SHM + isotype
```

---

## Modules

| Module | Tool | Input | Output |
|---|---|---|---|
| QC | FastQC + MultiQC | GEX + BCR FASTQ | HTML QC report |
| Alignment + GEX | Cell Ranger count | GEX FASTQ | BAM + gene matrix |
| TSS analysis | SCAFE | Cell Ranger BAM | tCRE matrix |
| Alternative TSS | CamoTSS | Cell Ranger BAM (filtered xf:i:25) | TC + CTSS |
| Enhancer calling | ReapTEC (STARsolo) | GEX FASTQ | btcEnh BED + DESeq2 |
| BCR repertoire | Cell Ranger VDJ | BCR FASTQ | Clonotypes |
| BCR repertoire | MiXCR | BCR FASTQ | SHM + isotype + clonal evolution |
| Clustering | Seurat 5.4 | Cell Ranger gene matrix | UMAP + cell type annotation |

---

## Requirements

- Nextflow ≥ 23.04
- Singularity (for HPC) or Docker
    - Cell Ranger (ref link: )
    - MiXCR (ref link: )
- ≥ 64GB RAM per sample
- ≥ 500GB disk space for full cohort

---

## Installation

```bash
git clone https://github.com/maitnnguyen/5prime-scrna
cd 5prime-scrna
```

Pre-pull Singularity images (run once on login node):

```bash
bash bin/pull_singularity_images.sh
```

---

## Samplesheet

Prepare a CSV file with one row per sample per lane:

```csv
sample,lane,condition,gex_r1,gex_r2,bcr_r1,bcr_r2
Control1,L001,Control,/path/Control1_S1_L001_R1_001.fastq.gz,/path/Control1_S1_L001_R2_001.fastq.gz,/path/Control1_BCR_S9_L001_R1_001.fastq.gz,/path/Control1_BCR_S9_L001_R2_001.fastq.gz
Control1,L002,Control,/path/Control1_S1_L002_R1_001.fastq.gz,/path/Control1_S1_L002_R2_001.fastq.gz,/path/Control1_BCR_S9_L002_R1_001.fastq.gz,/path/Control1_BCR_S9_L002_R2_001.fastq.gz
Patient1,L001,CVID,/path/Patient1_S5_L001_R1_001.fastq.gz,/path/Patient1_S5_L001_R2_001.fastq.gz,/path/Patient1_BCR_S13_L001_R1_001.fastq.gz,/path/Patient1_BCR_S13_L001_R2_001.fastq.gz
```

| Column | Description |
|---|---|
| `sample` | Sample ID (e.g. Control1, Patient1) |
| `lane` | Sequencing lane (L001 or L002) |
| `condition` | CVID or Control |
| `gex_r1` | Full path to GEX R1 FASTQ |
| `gex_r2` | Full path to GEX R2 FASTQ |
| `bcr_r1` | Full path to BCR R1 FASTQ |
| `bcr_r2` | Full path to BCR R2 FASTQ |

Generate samplesheet automatically:

```bash
python3 bin/generate_samplesheet.py
```

---

## Usage

### Run on SLURM cluster (Singularity)

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --outdir results/ \
    --genome hg38 \
    --cellranger_ref /path/to/refdata-gex-GRCh38-2020-A \
    --cellranger_vdj_ref /path/to/refdata-cellranger-vdj-GRCh38 \
    --scafe_genome hg38.gencode_v49 \
    --scafe_genome_dir /path/to/SCAFE_genome \
    --scafe_glm_model /path/to/hg38.gencode_v32/glm_model/SCAFE.glm.model.rds \
    --sif_dir /home/arkku/group/ics/tools/singularity \
    -profile singularity,slurm \
    -resume
```

### Run with Docker (local)

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --outdir results/ \
    --genome hg38 \
    --cellranger_ref /path/to/refdata-gex-GRCh38-2020-A \
    -profile docker \
    -resume
```

---

## Key Parameters

| Parameter | Default | Description |
|---|---|---|
| `--input` | required | Path to samplesheet CSV |
| `--outdir` | `results/` | Output directory |
| `--genome` | `hg38` | Reference genome |
| `--cellranger_ref` | required | Cell Ranger GEX reference folder |
| `--cellranger_vdj_ref` | required | Cell Ranger VDJ reference folder |
| `--scafe_genome` | `hg38.gencode_v49` | SCAFE genome name |
| `--scafe_genome_dir` | required | Path to SCAFE genome folder |
| `--scafe_glm_model` | required | Path to pre-trained GLM model (.rds) |
| `--sif_dir` | required | Path to folder containing pre-pulled .sif files |
| `--clip5p` | `39` | Bases to clip from R1 5' (Next GEM). **Change to 41 for GEM-X** |
| `--umi_len` | `12` | UMI length |
| `--bidir_gap` | `500` | Max gap (bp) between sense/antisense TSS pairs for ReapTEC |
| `--min_cpm` | `2` | log2CPM cutoff for TSS peak filtering |
| `--resolution` | `0.5` | Seurat clustering resolution |
| `--skip_scafe` | `false` | Skip SCAFE branch |
| `--skip_camotss` | `false` | Skip CamoTSS branch |
| `--skip_reaptec` | `false` | Skip ReapTEC branch |
| `--skip_vdj` | `false` | Skip Cell Ranger VDJ branch |
| `--skip_mixcr` | `false` | Skip MiXCR branch |

---

## GEM-X Chemistry

If using **10x GEM-X** chemistry (newer kit), adjust:

```bash
--clip5p 41 --umi_len 12
```

---

## Output Structure

```
results/
├── multiqc/                  ← Aggregated QC report (FastQC + Cell Ranger)
├── cellranger/
│   ├── {sample}/outs/        ← Cell Ranger count outputs
│   │   ├── possorted_genome_bam.bam
│   │   ├── filtered_feature_bc_matrix/
│   │   └── web_summary.html
│   └── vdj/{sample}/outs/    ← Cell Ranger VDJ outputs
├── mixcr/
│   └── {sample}/             ← MiXCR clonotypes + SHM + isotype
├── scafe/
│   └── {sample}/             ← tCRE BED + UMI count matrix
├── camotss/
│   └── {sample}/             ← TC + CTSS outputs
├── reaptec/
│   ├── {sample}/enhancers/   ← btcEnh BED files
│   └── pseudobulk/           ← DESeq2 CVID vs Control results
└── seurat/
    ├── {sample}_seurat.rds   ← Seurat objects
    └── {sample}_umap.pdf     ← UMAP plots
```

---

## Singularity Images

All containers are pre-pulled to avoid runtime downloads on compute nodes.
Images are stored in `--sif_dir` (default: `/home/arkku/group/ics/tools/singularity/`):

| File | Tool | Version |
|---|---|---|
| `scafe.sif` | [SCAFE](https://github.com/chung-lab/SCAFE) | latest |
| `fastqc.sif` | FastQC | 0.12.1 |
| `multiqc.sif` | MultiQC | 1.21 |
| `umi_tools.sif` | UMI-tools | 1.1.5 |
| `mixcr_latest.sif` | [MiXCR](https://mixcr.com/mixcr/getting-started/docker/) | 4.0.0 |
| `cellranger.sif` | [Cell Ranger](https://hub.docker.com/r/litd/docker-cellranger) | 8.0.1 |

**Seurat** (v5.4) is run via conda environment (`vscenv`) — no SIF required.
**samtools** and **STAR** are loaded via HPC modules.

Pull all images:
```bash
bash bin/pull_singularity_images.sh
```

---

## HPC Notes (arkku cluster)

```bash
# Request compute node
srun -p arkku --cpus-per-task=8 --mem=32G --time=08:00:00 --pty bash

# Activate conda environment
conda activate vscenv

# Run pipeline
nextflow run main.nf -profile singularity,slurm -resume
```

SCAFE requires custom genome binding — handled automatically via `--scafe_genome_dir`.

---

## Downstream Analysis

After running this pipeline, recommended next steps:

1. **TF Motif Enrichment** — chromVAR + JASPAR2024 on btcEnh BED
2. **GRN Reconstruction** — SCENIC+ using btcEnh + GEX
3. **BCR-TSS Integration** — link SCAFE tCRE activity to BCR isotype/SHM (via MiXCR)
4. **CVID Candidate Enhancers** — overlap DESeq2 results with GWAS loci
5. **Alternative TSS usage** — BRIE2 on CamoTSS output for differential TSS analysis

---

## Citation

If you use this pipeline, please cite:

- Oguchi A. et al. (2024). *Science* 385(6704). https://doi.org/10.1126/science.add8394
- Moody J. et al. (2022). SCAFE. *Bioinformatics*. https://doi.org/10.1093/bioinformatics/btac450
- Hou R. et al. (2023). CamoTSS. *Nat Commun* 14, 7240. https://doi.org/10.1038/s41467-023-42636-1
- Hao Y. et al. (2024). Seurat v5. *Nature Methods*. https://doi.org/10.1038/s41592-024-02353-z

---

## License

Apache-2.0 — see [LICENSE](LICENSE)
To run MiXCR, we obtained license for academic use from MiXCR (https://mixcr.com/mixcr/getting-started/milm/)
