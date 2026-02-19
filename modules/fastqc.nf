// ============================================================
//  modules/fastqc.nf
// ============================================================

process FASTQC {

    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}/fastqc/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip"),  emit: zip
    tuple val(meta), path("*.log"),  emit: log

    script:
    """
    fastqc \\
        --threads ${task.cpus} \\
        --outdir . \\
        ${reads} \\
        2> ${meta.id}.fastqc.log
    """
}
