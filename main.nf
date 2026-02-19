#!/usr/bin/env nextflow
// ============================================================
//  ReapTEC-B  |  B Cell Bidirectional Enhancer Pipeline  v2
//  ─────────────────────────────────────────────────────────
//  ReapTEC : Oguchi et al., Science 2024
//  SCAFE   : Moody & Kouno et al., Bioinformatics 2022
//  scTSS   : Fu & Li, GitHub 2024
//  ─────────────────────────────────────────────────────────
//  FASTQ → QC → STARsolo
//    ├──► ReapTEC (cap G → CTSS → btcEnh BED)
//    └──► SCAFE   (ML filter → tCRE BED)
//              ↓ cross-validate (bedtools intersect)
//         Consensus BED (TIER1 both / TIER2 ReapTEC / TIER3 SCAFE)
//              ↓
//         Seurat clustering · DESeq2 · scTSS · Outputs
// ============================================================

nextflow.enable.dsl = 2

// ── Import modules ───────────────────────────────────────────
include { FASTQC                  } from './modules/fastqc'
include { TRIMGALORE              } from './modules/trimgalore'
include { STAR_INDEX              } from './modules/star_index'
include { STAR_ALIGN              } from './modules/star_align'
include { UMITOOLS_DEDUP          } from './modules/umitools_dedup'
include { SOFTCLIP_G_FILTER       } from './modules/softclip_g'
include { CELL_BARCODE_CTSS       } from './modules/cell_barcode_ctss'
include { CTSS_COUNTS_BIGWIG      } from './modules/ctss_counts_bigwig'
include { BIDIR_ENHANCER_CALL     } from './modules/bidir_enhancer_call'
include { PROMOTER_FILTER         } from './modules/bidir_enhancer_call'
include { SCAFE_SOLO              } from './modules/scafe'
include { SCAFE_AGGREGATE         } from './modules/scafe'
include { SCAFE_COUNT_AGGREGATE   } from './modules/scafe'
include { CROSSVALIDATE_ENHANCERS } from './modules/crossvalidate'
include { SEURAT_CLUSTERING       } from './modules/seurat_clustering'
include { PSEUDOBULK_ENHANCER     } from './modules/seurat_clustering'
include { MULTIQC                 } from './modules/multiqc'

// ── Import subworkflows ──────────────────────────────────────
include { QC_SUBWORKFLOW    } from './subworkflows/qc'
include { REAPTEC_CORE      } from './subworkflows/reaptec_core'
include { SCAFE_WORKFLOW    } from './subworkflows/scafe_workflow'
include { ENHANCER_ANALYSIS } from './subworkflows/enhancer_analysis'

// ── Help ─────────────────────────────────────────────────────
def helpMessage() {
    log.info """
    ╔══════════════════════════════════════════════════════════════╗
    ║    ReapTEC-B v2  ·  B Cell Enhancer Pipeline                ║
    ║    ReapTEC + SCAFE + scTSS  ·  5' scRNA-seq                 ║
    ╚══════════════════════════════════════════════════════════════╝

    Usage:
        nextflow run main.nf --input samplesheet.csv \\
            --genome hg38 --scafe_genome hg38.gencode_v41 --outdir results/

    Required:
        --input             Samplesheet CSV (sample,fastq_R1,fastq_R2,condition)
        --genome            Genome build for STAR: hg38 (default) or mm10
        --scafe_genome      SCAFE genome name: hg38.gencode_v41
        --outdir            Output directory

    Reference files:
        --star_index        Pre-built STAR index (auto-built if absent)
        --gtf               Gene annotation GTF
        --fasta             Genome FASTA
        --whitelist         10x barcode whitelist (3M-february-2018.txt)
        --fantom_bed        FANTOM5 promoter/enhancer BED (auto-downloaded)

    ReapTEC parameters:
        --clip5p            Bases to clip from R1 5' end: 39 (Next GEM) | 41 (GEM-X)
        --umi_len           UMI length: 12
        --bidir_gap         Max gap (bp) for bidirectional pairing: 500
        --min_cpm           log2CPM cutoff for robust TSS: 2
        --min_mapq          Minimum MAPQ: 255

    SCAFE parameters:
        --scafe_genome      SCAFE genome name (must match downloaded reference)
        --scafe_atac_bw     (Optional) ATAC bigwig for SCAFE logistic model training
        --scafe_glm_model   (Optional) Pre-built GLM model (skips ATAC training)
        --skip_scafe        Skip SCAFE entirely [false]
        --skip_crossval     Skip cross-validation step [false]

    Analysis parameters:
        --min_cells         Min cells per feature (Seurat): 3
        --min_features      Min features per cell: 200
        --max_mt_pct        Max mitochondrial % cutoff: 20
        --resolution        Seurat clustering resolution: 0.5
        --n_pcs             Number of PCs: 30
        --skip_vdj          Skip BCR V(D)J extraction [false]

    Compute:
        --max_cpus / --max_memory / --max_time

    Profiles:
        -profile docker / singularity / slurm / test

    Examples:
        # Full pipeline (ReapTEC + SCAFE + cross-validate)
        nextflow run main.nf --input samplesheet.csv \\
            --genome hg38 --scafe_genome hg38.gencode_v41 \\
            --outdir results/ -profile docker

        # ReapTEC only (skip SCAFE)
        nextflow run main.nf --input samplesheet.csv \\
            --genome hg38 --skip_scafe true \\
            --outdir results/ -profile docker
    """.stripIndent()
}

if (params.help) { helpMessage(); exit 0 }

// ── Validate ─────────────────────────────────────────────────
if (!params.input)  { error "ERROR: --input samplesheet is required" }
if (!params.outdir) { error "ERROR: --outdir is required" }
if (!params.skip_scafe && !params.scafe_genome) {
    error "ERROR: --scafe_genome is required for SCAFE (or set --skip_scafe true)"
}

// ── Main workflow ─────────────────────────────────────────────
workflow {

    log.info """
    ╔══════════════════════════════════════════════════════════════╗
    ║  ReapTEC-B v2  |  Starting Run                              ║
    ╠══════════════════════════════════════════════════════════════╣
    ║  Input          : ${params.input}
    ║  STAR genome    : ${params.genome}
    ║  SCAFE genome   : ${params.scafe_genome ?: 'SKIPPED (--skip_scafe)'}
    ║  Output         : ${params.outdir}
    ║  ReapTEC clip5p : ${params.clip5p} bp
    ║  Bidir gap      : ${params.bidir_gap} bp
    ║  Skip SCAFE     : ${params.skip_scafe}
    ║  Skip CrossVal  : ${params.skip_crossval}
    ╚══════════════════════════════════════════════════════════════╝
    """.stripIndent()

    // ── 1. Parse samplesheet ──────────────────────────────────
    ch_samplesheet = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def meta = [
                id        : row.sample,
                condition : row.condition ?: 'unknown',
                single_end: false
            ]
            [ meta, [ file(row.fastq_R1, checkIfExists: true),
                      file(row.fastq_R2, checkIfExists: true) ] ]
        }

    // ── 2. Reference channels ─────────────────────────────────
    ch_fasta     = params.fasta     ? Channel.value(file(params.fasta))     : Channel.empty()
    ch_gtf       = params.gtf       ? Channel.value(file(params.gtf))       : Channel.empty()
    ch_whitelist = Channel.value(file(params.whitelist))
    ch_fantom    = params.fantom_bed ? Channel.value(file(params.fantom_bed)) : Channel.empty()

    if (!params.star_index) {
        STAR_INDEX(ch_fasta, ch_gtf)
        ch_star_index = STAR_INDEX.out.index
    } else {
        ch_star_index = Channel.value(file(params.star_index))
    }

    // ── 3. QC ─────────────────────────────────────────────────
    QC_SUBWORKFLOW(ch_samplesheet)
    ch_trimmed = QC_SUBWORKFLOW.out.reads
    ch_qc_logs = QC_SUBWORKFLOW.out.logs

    // ── 4. STARsolo alignment (shared BAM → both branches) ────
    REAPTEC_CORE(ch_trimmed, ch_star_index, ch_gtf, ch_whitelist)
    ch_bam      = REAPTEC_CORE.out.bam
    ch_ctss     = REAPTEC_CORE.out.ctss_bed
    ch_barcodes = REAPTEC_CORE.out.barcodes
    ch_gex      = REAPTEC_CORE.out.gene_matrix

    // ── 5a. ReapTEC branch ────────────────────────────────────
    // Cap G filter → per-cell CTSS → bidirectional enhancer call
    BIDIR_ENHANCER_CALL(ch_ctss, ch_fantom)
    PROMOTER_FILTER(BIDIR_ENHANCER_CALL.out.raw_enhancers, ch_fantom)
    ch_reaptec_enhancers = PROMOTER_FILTER.out.enhancers

    // ── 5b. SCAFE branch (parallel, same BAM) ─────────────────
    // ML logistic regression → tCRE BED → UMI count matrix
    if (!params.skip_scafe) {

        SCAFE_WORKFLOW(ch_bam, ch_bam, ch_barcodes)
        ch_scafe_tcre     = SCAFE_WORKFLOW.out.tCRE_bed
        ch_scafe_logs     = SCAFE_WORKFLOW.out.logs
        ch_scafe_matrices = SCAFE_WORKFLOW.out.count_matrices

        // ── 6. Cross-validate ReapTEC × SCAFE ─────────────────
        // TIER 1: overlap both tools → highest confidence
        // TIER 2: ReapTEC only (bidirectional, not in SCAFE)
        // TIER 3: SCAFE only (ML validated, not bidirectional)
        if (!params.skip_crossval) {
            CROSSVALIDATE_ENHANCERS(
                ch_reaptec_enhancers,
                ch_reaptec_enhancers.combine(ch_scafe_tcre.first())
                    .map { meta, reaptec, scafe -> [ meta, scafe ] }
            )
            ch_consensus_bed = CROSSVALIDATE_ENHANCERS.out.consensus
            ch_crossval_logs = CROSSVALIDATE_ENHANCERS.out.log
        } else {
            ch_consensus_bed = ch_reaptec_enhancers
            ch_crossval_logs = Channel.empty()
        }

    } else {
        ch_consensus_bed  = ch_reaptec_enhancers
        ch_scafe_logs     = Channel.empty()
        ch_scafe_matrices = Channel.empty()
        ch_crossval_logs  = Channel.empty()
        log.warn "SCAFE branch skipped. Using ReapTEC btcEnhs only."
    }

    // ── 7. Downstream analysis ────────────────────────────────
    // Seurat B cell clustering + DESeq2 pseudobulk
    ENHANCER_ANALYSIS(ch_ctss, ch_fantom, ch_bam, ch_gex, ch_barcodes)

    // ── 8. MultiQC ────────────────────────────────────────────
    ch_all_logs = ch_qc_logs
        .mix(REAPTEC_CORE.out.logs)
        .mix(ENHANCER_ANALYSIS.out.logs)
        .mix(ch_scafe_logs)
        .mix(ch_crossval_logs)
        .collect()

    MULTIQC(ch_all_logs)

    // ── 9. Emit final outputs ─────────────────────────────────
    emit:
    consensus_enhancers = ch_consensus_bed
    reaptec_btcEnhs     = ch_reaptec_enhancers
    gex_matrix          = ch_gex
    diff_enhancers      = ENHANCER_ANALYSIS.out.diff_enhancers
    multiqc_report      = MULTIQC.out.report
}

workflow.onComplete {
    log.info """
    ╔══════════════════════════════════════════════════════════════╗
    ║  Pipeline Complete!                                          ║
    ║  Status   : ${workflow.success ? 'SUCCESS ✓' : 'FAILED ✗'}
    ║  Duration : ${workflow.duration}
    ╠══════════════════════════════════════════════════════════════╣
    ║  results/crossvalidate/  ← Consensus BED (TIER 1/2/3)       ║
    ║  results/enhancers/      ← ReapTEC btcEnh BED               ║
    ║  results/scafe/          ← SCAFE tCRE BED + UMI matrix      ║
    ║  results/seurat/         ← B cell clusters + UMAP           ║
    ║  results/pseudobulk/     ← DESeq2 CVID vs HC                ║
    ║  results/multiqc/        ← QC report                        ║
    ╚══════════════════════════════════════════════════════════════╝
    """.stripIndent()
}
