#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ============================================================
//  main.nf  |  ReapTEC-B Pipeline
//  5' scRNA-seq → CellRanger → STARsolo/ReapTEC → Enhancers
//
//  Usage:
//    nextflow run main.nf                           # hg38, full run
//    nextflow run main.nf --genome t2t              # T2T genome
//    nextflow run main.nf --skip_cellranger         # reuse existing CellRanger output
//    nextflow run main.nf --skip_fastqc             # skip FastQC/MultiQC
//    nextflow run main.nf -profile singularity      # use containers
// ============================================================

include { CELLRANGER_WORKFLOW } from './subworkflows/cellranger'
include { REAPTEC_CORE        } from './subworkflows/reaptec_core'
include { EXTRACT_WHITELIST   } from './modules/extract_whitelist'
include { FASTQC              } from './modules/fastqc'
include { MULTIQC             } from './modules/multiqc'

// ── Validate inputs ──────────────────────────────────────────
if (!params.input) {
    error "Please provide --input /path/to/samplesheet.csv"
}
if (!params.genomes.containsKey(params.genome)) {
    error "Genome '${params.genome}' not found. Valid options: ${params.genomes.keySet().join(', ')}"
}

// ── Resolve genome references once ───────────────────────────
// All paths resolved here and passed explicitly to subworkflows.
// Never scatter params.* lookups into modules — hard to trace.
def g = params.genomes[params.genome]

def star_index   = file(g.star_index,  checkIfExists: true)
def gtf          = file(g.gtf,         checkIfExists: true)
def gex_ref      = file(g.gex_ref,     checkIfExists: true)
def vdj_ref      = file(g.vdj_ref,     checkIfExists: true)
def mask_bed     = file(g.mask_bed,    checkIfExists: true)
def chrom_sizes  = file(g.chrom_sizes, checkIfExists: true)

// FANTOM references optional — not available for T2T
def fantom_promo = g.fantom_promo ? file(g.fantom_promo) : []
def fantom_enhan = g.fantom_enhan ? file(g.fantom_enhan) : []

// ── Log run info ─────────────────────────────────────────────
log.info """
    ╔══════════════════════════════════════════╗
    ║          ReapTEC-B Pipeline              ║
    ╚══════════════════════════════════════════╝
    input           : ${params.input}
    genome          : ${params.genome}
    outdir          : ${params.outdir}
    skip_cellranger : ${params.skip_cellranger}
    skip_fastqc     : ${params.skip_fastqc}
    solo_strand     : ${params.solo_strand}
    umi_len         : ${params.umi_len}
    clip5p          : ${params.clip5p}
    """.stripIndent()

// ============================================================
workflow {

    // ── Parse samplesheet ─────────────────────────────────────
    // Required columns : sample, fastq_dir
    // Optional columns : lane (default L001), cellranger_dir
    //
    // cellranger_dir is only needed when --skip_cellranger is set.
    // It should point to an existing CellRanger count output directory,
    // e.g. /path/to/results/hg38/cellranger/SAMPLE1
    //
    // Example samplesheet:
    //   sample,fastq_dir,lane,cellranger_dir
    //   SAMPLE1,/data/fastqs,L001,
    //   SAMPLE2,/data/fastqs,L001,/results/hg38/cellranger/SAMPLE2

    Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            if (!row.sample)    error "Samplesheet missing 'sample' column"
            if (!row.fastq_dir) error "Samplesheet missing 'fastq_dir' column for sample ${row.sample}"

            def meta           = [ id: row.sample ]   // ← remove lane
            def fastq_dir      = file(row.fastq_dir, checkIfExists: true)
            def cellranger_dir = row.cellranger_dir?.trim() ?: ''

            return [ meta, fastq_dir, cellranger_dir ]
        }
        .set { ch_samples }

    // ── Reference channels ────────────────────────────────────
    ch_star_index   = Channel.value(star_index)
    ch_gtf          = Channel.value(gtf)
    ch_gex_ref      = Channel.value(gex_ref)
    ch_vdj_ref      = Channel.value(vdj_ref)
    ch_mask_bed     = Channel.value(mask_bed)
    ch_chrom_sizes  = Channel.value(chrom_sizes)
    ch_fantom_promo = Channel.value(fantom_promo)
    ch_fantom_enhan = Channel.value(fantom_enhan)

    // ── FastQC (optional, independent of CellRanger) ─────────
    // Runs per-file per-lane on both GEX and BCR FASTQs.
    // Per-lane is intentional: lets you spot lane-specific QC issues
    // that would be hidden if you merged before QC.
    if (!params.skip_fastqc) {
        ch_fastqc_input = ch_samples
            .flatMap { meta, fastq_dir, cellranger_dir ->
                // Find all GEX and BCR FASTQs for this sample
                def gex_files = file("${fastq_dir}/${meta.id}_S*_L*_R*_001.fastq.gz")
                def bcr_files = file("${fastq_dir}/${meta.id}_BCR_S*_L*_R*_001.fastq.gz")
                (gex_files + bcr_files).collect { f -> [ meta, f ] }
            }

        FASTQC(ch_fastqc_input)

        // Collect all FastQC zips + CellRanger metrics for MultiQC
        // CellRanger metrics added below after CellRanger block
        ch_multiqc_input = FASTQC.out.zip
            .map { meta, zip -> zip }
            .collect()

        MULTIQC(ch_multiqc_input)
    }

    // ── CellRanger: run or skip ───────────────────────────────
    if (!params.skip_cellranger) {

        // Both GEX and BCR use the same fastq_dir —
        // CellRanger distinguishes them via --sample flag:
        //   COUNT: --sample=SAMPLE1       matches SAMPLE1_S*_R*
        //   VDJ:   --sample=SAMPLE1_BCR   matches SAMPLE1_BCR_S*_R*
        ch_gex_reads = ch_samples.map { meta, fastq_dir, cellranger_dir -> [ meta, fastq_dir ] }
        ch_bcr_reads = ch_samples.map { meta, fastq_dir, cellranger_dir -> [ meta, fastq_dir ] }

        CELLRANGER_WORKFLOW(
            ch_gex_reads,
            ch_bcr_reads,
            ch_gex_ref,
            ch_vdj_ref
        )

        ch_matrix = CELLRANGER_WORKFLOW.out.matrix

    } else {

        // Use existing CellRanger output — cellranger_dir must be set
        // in samplesheet for every sample when this flag is used
        ch_matrix = ch_samples
        .map { meta, fastq_dir, cellranger_dir ->
            if (!cellranger_dir) {
                error "Sample ${meta.id}: --skip_cellranger is set but 'cellranger_dir' is missing in samplesheet"
            }
            def matrix = file("${cellranger_dir}/outs/filtered_feature_bc_matrix", checkIfExists: true)
            return [ meta, matrix ]
        }
        .unique { meta, matrix -> meta.id }  // deduplicate — one matrix per sample
    }

    // ── Extract per-sample whitelist ──────────────────────────
    // Always runs — works the same whether matrix came from a fresh
    // CellRanger run or from an existing output directory
    EXTRACT_WHITELIST(ch_matrix)

    // ── Join reads + whitelist for STAR ───────────────────────
    // Join by meta.id to guarantee correct whitelist per sample
    ch_star_input = ch_samples
        .map    { meta, fastq_dir, cellranger_dir -> [ meta, fastq_dir ] }
        .unique { meta, fastq_dir -> meta.id }   // ← distict the path due to multiple lanes per sample
        .join   ( EXTRACT_WHITELIST.out.whitelist, by: 0 )
        // → [ meta, fastq_dir, whitelist ]
        // Remap to reads tuple expected by STAR_ALIGN:
        // STAR needs actual R1/R2 files, not just the directory
        .flatMap { meta, fastq_dir, whitelist ->
            def r1_files = file("${fastq_dir}/${meta.id}_S*_L*_R1_001.fastq.gz").sort()
            def r2_files = file("${fastq_dir}/${meta.id}_S*_L*_R2_001.fastq.gz").sort()
            if (!r1_files) error "No R1 FASTQs found for sample ${meta.id} in ${fastq_dir}"
            if (!r2_files) error "No R2 FASTQs found for sample ${meta.id} in ${fastq_dir}"
            // Return: [ meta, r1_files, r2_files, whitelist ] — 4 elements
            // STAR accepts comma-separated / multiple files for multi-lane
            return [ [ meta, r1_files, r2_files, whitelist ] ]
        }

    // ── ReapTEC core ──────────────────────────────────────────
    REAPTEC_CORE(
        ch_star_input,
        ch_star_index,
        ch_gtf,
        ch_chrom_sizes,
        ch_fantom_promo,
        ch_fantom_enhan,
        ch_mask_bed
    )

    // ── Completion messages ───────────────────────────────────
    REAPTEC_CORE.out.ctss_bed.view { meta, bed ->
        "[ReapTEC] ${meta.id} (${params.genome}): CTSS complete → ${bed}"
    }
    REAPTEC_CORE.out.enhancers.view { bed ->
        "[ReapTEC] Enhancers called → ${bed}"
    }
}
