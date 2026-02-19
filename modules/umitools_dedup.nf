// ============================================================
//  modules/umitools_dedup.nf
//  UMI-based PCR duplicate removal
//  Must run AFTER STARsolo so CB/UB tags are present in BAM
// ============================================================

process UMITOOLS_DEDUP {

    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/umitools/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(bam)
    tuple val(meta), path(bai)

    output:
    tuple val(meta), path("${meta.id}.dedup.bam"),     emit: bam
    tuple val(meta), path("${meta.id}.dedup.bam.bai"), emit: bai
    tuple val(meta), path("${meta.id}.dedup.log"),     emit: log

    script:
    """
    umi_tools dedup \\
        --stdin=${bam} \\
        --stdout=${meta.id}.dedup.bam \\
        --per-cell \\
        --cell-tag=CB \\
        --umi-tag=UB \\
        --extract-umi-method=tag \\
        --method=unique \\
        --log=${meta.id}.dedup.log

    samtools index -@ ${task.cpus} ${meta.id}.dedup.bam
    """
}
