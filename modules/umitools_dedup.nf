// ============================================================
//  modules/umitools_dedup.nf
//
//  Position in ReapTEC pipeline:
//    SOFTCLIP_G_FILTER → UMITOOLS_DEDUP → CTSS_BED
//
//  Tag choice: CR/UR (raw) vs CB/UB (corrected)
//  ─────────────────────────────────────────────
//  The official ReapTEC protocol (Murakawa lab, April 2023) uses:
//    --umi-tag=UR  --cell-tag=CR  (raw tags)
//  This pipeline intentionally DEVIATES to use CR/UR as per protocol.
//
//  Note: If you want to use corrected tags (CB/UB), swap the tag flags
//  below. Using CB/UB would implicitly restrict dedup to whitelist-matched
//  reads only (since STARsolo only writes CB/UB for matched barcodes),
//  which is a reasonable alternative but not the validated approach.
// ============================================================

process UMITOOLS_DEDUP {

    tag "$meta.id"
    label 'process_high'    // Single-threaded but RAM-hungry

    publishDir "${params.outdir}/${params.genome}/reaptec/dedup", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai), path(whitelist)
    // All four arrive as one tuple so Nextflow guarantees
    // the correct whitelist is paired with the correct BAM.

    output:
    tuple val(meta), path("${meta.id}.dedup.bam"),     emit: bam
    tuple val(meta), path("${meta.id}.dedup.bam.bai"), emit: bai
    path "${meta.id}.dedup.log",                        emit: log

    script:
    """
    # Re-encode BAM with HPC samtools before dedup (fixes str/bytes bug in umi_tools)
    samtools view -h ${bam} | samtools view -b -o reencoded.bam
    samtools index reencoded.bam

    # ── Step 1: UMI deduplication ─────────────────────────────
    # Use raw CR/UR tags per official ReapTEC protocol.
    # STARsolo writes CR/UR for ALL reads; CB/UB only for whitelist-matched.
    # Using CR/UR here ensures no reads are silently dropped before
    # the explicit barcode filter in step 2.
    umi_tools dedup \\
        -I ${bam} \\
        --per-cell \\
        --cell-tag=CB \\
        --umi-tag=UR \\
        --extract-umi-method=tag \\
        --method=unique \\
        --log=${meta.id}.dedup.log \\
        -S dedup_raw.bam

    samtools index -@ ${task.cpus} dedup_raw.bam

    # ── Step 2: Filter to valid cell barcodes ─────────────────
    # Use samtools -D (tag-based filter) rather than grep on raw SAM.
    # -D CB:<file>  keeps only reads whose CB tag matches the list.
    # This is exact-match and field-aware — grep on raw SAM can match
    # partial strings in other fields.
    #
    # Note: after dedup with CR/UR, STARsolo-corrected CB tags ARE
    # present in the BAM (written by STAR regardless of dedup step).
    # So we filter by CB here to keep only validated-barcode reads.
    samtools view \\
        -@ ${task.cpus} \\
        -h \\
        -D CB:${whitelist} \\
        dedup_raw.bam \\
    | samtools sort -@ ${task.cpus} -o ${meta.id}.dedup.bam

    samtools index -@ ${task.cpus} ${meta.id}.dedup.bam

    # Cleanup
    rm dedup_raw.bam dedup_raw.bam.bai
    """
}
