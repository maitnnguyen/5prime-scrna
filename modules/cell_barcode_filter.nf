// ============================================================
//  modules/filter_barcodes.nf
//  Filter deduplicated BAM to valid cell barcodes only.
//
//  Uses samtools -D CB to keep only reads whose CB tag
//  matches the per-sample CellRanger whitelist.
//  This is exact tag-value matching — safer than grep on raw SAM.
//
//  No container — uses HPC samtools module.
// ============================================================

process CELL_BARCODE_FILTER {

    tag "$meta.id"
    label 'process_medium'

    publishDir "${params.outdir}/${params.genome}/reaptec/dedup", mode: 'copy'

    input:
    tuple val(meta), path(bam)
    tuple val(meta), path(whitelist)

    output:
    tuple val(meta), path("${meta.id}.dedup.filtered.bam"),     emit: bam
    tuple val(meta), path("${meta.id}.dedup.filtered.bam.bai"), emit: bai

    script:
    """
    # Filter to valid cell barcodes using CB tag (corrected barcode)
    # -D CB:<file> keeps only reads whose CB tag matches the whitelist
    # This is field-aware exact matching — safer than grep
    samtools view \\
        -@ ${task.cpus} \\
        -h \\
        -D CB:${whitelist} \\
        ${bam} \\
    | samtools sort -@ ${task.cpus} -o ${meta.id}.dedup.filtered.bam

    samtools index -@ ${task.cpus} ${meta.id}.dedup.filtered.bam

    # Report final cell count
    NCELLS=\$(samtools view -@ ${task.cpus} ${meta.id}.dedup.filtered.bam | \\
        grep -o 'CB:Z:[^ ]*' | sort -u | wc -l)
    echo "[FILTER_BARCODES] ${meta.id}: \$NCELLS unique cell barcodes retained"
    """
}
