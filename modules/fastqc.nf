// ============================================================
//  modules/fastqc.nf
//  FastQC — Per-sample read quality control
//
//  Notes:
//  - Runs on both GEX and BCR FASTQ files
//  - HTML/zip reports only — no FASTQ files saved (space saving)
//  - MultiQC reads zip files directly, no need to extract
//  - --noextract keeps zip compressed (MultiQC compatible)
//  - --nogroup gives per-base quality, not grouped (more detail)
// ============================================================

process FASTQC {

    tag "${meta.id} ${meta.library} ${meta.lane}"
    label 'process_low'

    // Only publish QC reports — not FASTQ files
    publishDir "${params.outdir}/fastqc/${meta.id}", mode: 'copy',
        saveAs: { filename -> filename.endsWith('.log') ? null : filename }

    input:
    tuple val(meta), path(reads)   // meta has: id, library (GEX/BCR), lane

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip"),  emit: zip

    script:
    """
    fastqc \\
        --threads ${task.cpus} \\
        --outdir . \\
        --noextract \\
        --nogroup \\
        ${reads} \\
        2> ${meta.id}_${meta.library}_${meta.lane}.fastqc.log
    """
}
