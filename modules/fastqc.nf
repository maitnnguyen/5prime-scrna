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
    tag "${sample_id}"
    label 'process_low'
    publishDir "${params.outdir}/fastqc", mode: 'copy', pattern: "*.html"

    container params.fastqc_container

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.html"), emit: html
    tuple val(sample_id), path("*.zip"),  emit: zip
    path  "*.zip",                        emit: reports

    script:
    """
    fastqc \\
        --threads ${task.cpus} \\
        --quiet \\
        --outdir . \\
        ${reads}
    """
}