// ============================================================
//  modules/umitools_dedup.nf
//  UMI deduplication using umi_tools.
//
//  Position in pipeline:
//    PREPARE_DEDUP_BAM → UMITOOLS_DEDUP → FILTER_BARCODES
//
//  Uses umi_tools Singularity container — no samtools inside.
//  samtools steps are handled by PREPARE_DEDUP_BAM and FILTER_BARCODES.
//
//  Tag choice: CR/UR (raw tags) per official ReapTEC protocol.
//  PREPARE_DEDUP_BAM ensures all UR tags are uniform length
//  before this step to avoid AssertionError in umi_tools.
// ============================================================

process UMITOOLS_DEDUP {

    tag "$meta.id"
    label 'process_high'   // Single-threaded but RAM-hungry

    publishDir "${params.outdir}/${params.genome}/reaptec/dedup", mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith('.log')) return filename
            else null   // BAM published after barcode filtering in FILTER_BARCODES
        }

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.dedup.bam"), emit: bam
    path "${meta.id}.dedup.log",                   emit: log

    script:
    """
    umi_tools dedup \\
        -I ${bam} \\
        --per-cell \\
        --extract-umi-method=tag \\
        --cell-tag=CR \\
        --umi-tag=UR \\
        --log=${meta.id}.dedup.log \\
        -S ${meta.id}.dedup.bam
    """
}
