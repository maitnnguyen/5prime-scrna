// ============================================================
//  modules/trimgalore.nf
//  Adapter trimming for 5' scRNA-seq paired-end reads
// ============================================================

process TRIMGALORE {

    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/trimgalore/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_val_{1,2}.fq.gz"),         emit: reads
    tuple val(meta), path("*_trimming_report.txt"),      emit: log
    tuple val(meta), path("*.html"),                     emit: html, optional: true

    script:
    def r1 = reads[0]
    def r2 = reads[1]
    """
    trim_galore \\
        --paired \\
        --gzip \\
        --cores ${task.cpus} \\
        --fastqc \\
        ${r1} ${r2}
    """
}
