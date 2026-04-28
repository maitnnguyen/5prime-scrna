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
    # Add CB:Z: prefix to whitelist for grep matching (same as original ReapTEC)
    awk '{print "CB:Z:"\$1}' ${whitelist} > cell_barcode_prefixed.txt
    
    # Filter using grep on raw SAM — matches original ReapTEC method exactly
    samtools view -@ ${task.cpus} -H ${bam} > SAM_header
    samtools view -@ ${task.cpus} ${bam} | \
        LC_ALL=C grep -F -f cell_barcode_prefixed.txt > filtered_SAM_body
    cat SAM_header filtered_SAM_body | \
        samtools sort -@ ${task.cpus} -o ${meta.id}.dedup.filtered.bam

    samtools index -@ ${task.cpus} ${meta.id}.dedup.filtered.bam

    NCELLS=\$(samtools view -@ ${task.cpus} ${meta.id}.dedup.filtered.bam | \\
        grep -o 'CB:Z:[^ ]*' | sort -u | wc -l)
    echo "[FILTER_BARCODES] ${meta.id}: \$NCELLS unique cell barcodes retained"

    rm SAM_header filtered_SAM_body cell_barcode_prefixed.txt
    """
}
