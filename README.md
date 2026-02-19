# ============================================================
#  ReapTEC-B Pipeline | README
#  B Cell Bidirectional Enhancer Analysis from 5' scRNA-seq
# ============================================================

## Overview

This Nextflow pipeline implements the **ReapTEC** approach
(Oguchi et al., *Science* 2024) for identifying bidirectionally
transcribed candidate enhancers (btcEnhs) from 5' scRNA-seq data,
adapted specifically for **B cell / CVID research**.

**Reference:** Oguchi A. et al., Science 385(6704), 2024.
doi:10.1126/science.add8394

**Original ReapTEC:** github.com/MurakawaLab/ReapTEC
**nf-core ReapTEC:** github.com/paganilab/nf-core-reaptec

---

## Pipeline Overview

```
FASTQ (5' scRNA-seq)
    │
    ▼
[QC] FastQC → Trim Galore → MultiQC
    │
    ▼
[ALIGN] STARsolo (5' strand-aware, soloStrand=Forward)
    │   clip5pNbases=39 | CB:16bp | UMI:12bp
    ▼
[DEDUP] UMI-tools deduplication
    │
    ▼
[REAPTEC] SoftclipG Filter ──── KEY STEP
    │   Extract reads with 5' cap G signature
    │   Sense:    CIGAR starts with 39S + base G at pos 40
    │   Antisense: CIGAR ends with 39S + base C at 3' end
    ▼
[CTSS] Per-cell TSS BED generation
    │   chr  TSS_pos  cell_barcode  strand
    ▼
[BIGWIG] CPM-normalized BigWig (IGV/UCSC visualization)
    │
    ├──► [GEX] STARsolo Gene x Cell matrix → Seurat
    │         B cell clustering + subtype annotation
    │
    ▼
[ENHANCERS] Bidirectional Enhancer Calling
    │   1. Pseudo-bulk CTSS aggregation
    │   2. Peak clustering (merge within 10bp)
    │   3. Sense + antisense peak separation
    │   4. Bidirectional pairs (gap ≤ 500bp)
    │   5. Promoter masking (FANTOM5 + gene TSS ±2kb)
    │   6. CPM filter (log2CPM ≥ 2)
    ▼
[PSEUDOBULK] Per-cluster enhancer count matrix
    │   DESeq2: CVID vs Control differential analysis
    ▼
[OUTPUT] 
    ├── results/enhancers/   ← btcEnh BED files
    ├── results/seurat/      ← Seurat RDS + UMAP
    ├── results/pseudobulk/  ← Count matrix + DESeq2
    ├── results/bigwig/      ← IGV tracks
    └── results/multiqc/     ← QC report
```

---

## Quick Start

### 1. Install Nextflow
```bash
curl -s https://get.nextflow.io | bash
mv nextflow ~/bin/
```

### 2. Prepare samplesheet (samplesheet.csv)
```
sample,fastq_R1,fastq_R2,condition
CVID_patient1,/data/P1_R1.fastq.gz,/data/P1_R2.fastq.gz,CVID
CVID_patient2,/data/P2_R1.fastq.gz,/data/P2_R2.fastq.gz,CVID
HC_control1,/data/HC1_R1.fastq.gz,/data/HC1_R2.fastq.gz,Control
HC_control2,/data/HC2_R1.fastq.gz,/data/HC2_R2.fastq.gz,Control
```

### 3. Run pipeline (Docker)
```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --genome hg38 \
    --outdir results/ \
    -profile docker
```

### 4. Run on SLURM cluster (Singularity)
```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --genome hg38 \
    --outdir results/ \
    --max_cpus 32 \
    --max_memory 256.GB \
    -profile singularity,slurm \
    -resume
```

---

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--clip5p` | 39 | Bases to clip from R1 5' (Next GEM). **Change to 41 for GEM-X** |
| `--umi_len` | 12 | UMI length (Next GEM / GEM-X) |
| `--bidir_gap` | 500 | Max gap (bp) between sense/antisense TSS pairs |
| `--min_cpm` | 2 | log2CPM cutoff for robust TSS peaks |
| `--resolution` | 0.5 | Seurat clustering resolution |
| `--min_mapq` | 255 | Min MAPQ (255 = STAR uniquely mapped) |

---

## GEM-X Chemistry Adjustment

If using **10x GEM-X** chemistry (newer kit), change:
```bash
--clip5p 41 \
--umi_len 12
```
And the softclipG script automatically adjusts position indices.

---

## Output Structure

```
results/
├── pipeline_info/
│   ├── execution_report.html
│   ├── execution_timeline.html
│   └── pipeline_dag.svg
├── fastqc/          ← Raw read QC
├── trimgalore/      ← Trimming reports
├── star/            ← Aligned BAMs + STARsolo output
├── softclip_g/      ← Cap-signature filtered BAMs
├── ctss/            ← Per-cell CTSS BED files
├── bigwig/          ← Normalized tracks (sense + antisense)
├── enhancers/       ← btcEnh BED files per sample
│   ├── *.sense_peaks.bed
│   ├── *.antisense_peaks.bed
│   ├── *.bidir_pairs.bed
│   └── *.enhancers.bed      ← FINAL enhancers
├── seurat/          ← Seurat objects + UMAP plots
├── pseudobulk/      ← Count matrices + DESeq2 results
└── multiqc/         ← Aggregated QC report
```

---

## Downstream Analysis (Post-pipeline)

After running this pipeline, recommended next steps:

1. **TF Motif Enrichment** — chromVAR + JASPAR2024 on btcEnh BED
2. **GRN Reconstruction** — SCENIC+ using btcEnh + GEX
3. **BCR Integration** — Link btcEnh activity to isotype (Cell Ranger VDJ)
4. **CVID Candidate Enhancers** — Overlap DESeq2 results with GWAS loci
5. **Super-Enhancers** — ROSE on differentially active enhancers

---

## Citation

If you use this pipeline, please cite:

> Oguchi A. et al. (2024). "Single-cell analysis of human CD4+ T cells
> identifies enhancers associated with autoimmune diseases."
> *Science* 385(6704). https://doi.org/10.1126/science.add8394

> Andersson R. et al. (2014). "An atlas of active enhancers across human
> cell types and tissues." *Nature* 507, 455–461.
