// ============================================================
//  subworkflows/cellranger_multi.nf
//  CellRanger COUNT + VDJ only.
//
//  FastQC and MultiQC are intentionally NOT here —
//  they run independently in main.nf with --skip_fastqc flag.
// ============================================================

include { CELLRANGER_COUNT } from '../modules/cellranger_count'
include { CELLRANGER_VDJ   } from '../modules/cellranger_vdj'

workflow CELLRANGER_WORKFLOW {

    take:
    ch_gex_reads  // [ meta, fastq_dir ] — GEX samples
    ch_bcr_reads  // [ meta, fastq_dir ] — BCR samples (filtered to those with BCR FASTQs)
    ch_gex_ref    // CellRanger GEX reference (value channel)
    ch_vdj_ref    // CellRanger VDJ reference (value channel)

    main:

    // ── CellRanger count (GEX) ────────────────────────────────
    // --sample=${meta.id} matches SAMPLE1_S*_R* files only
    // BCR files (SAMPLE1_BCR_*) are automatically excluded
    CELLRANGER_COUNT(ch_gex_reads, ch_gex_ref)

    // ── CellRanger VDJ (BCR) ─────────────────────────────────
    // --sample=${meta.id}_BCR matches SAMPLE1_BCR_S*_R* files only
    // Only runs for samples that have BCR FASTQs (filtered in main.nf)
    CELLRANGER_VDJ(ch_bcr_reads, ch_vdj_ref)

    emit:
    bam          = CELLRANGER_COUNT.out.bam
    bai          = CELLRANGER_COUNT.out.bai
    matrix       = CELLRANGER_COUNT.out.matrix      // → EXTRACT_WHITELIST
    barcodes     = CELLRANGER_COUNT.out.barcodes
    summary      = CELLRANGER_COUNT.out.summary
    metrics_gex  = CELLRANGER_COUNT.out.metrics
    clonotypes   = CELLRANGER_VDJ.out.clonotypes
    airr         = CELLRANGER_VDJ.out.airr
    metrics_vdj  = CELLRANGER_VDJ.out.metrics
}
