// ============================================================
//  modules/cell_barcode_ctss.nf
//  ReapTEC Step 4: Per-cell CTSS BED generation
//  Adapted from: STARsolo_Cell_barcode_CTSS_bed_20221123.sh
// ============================================================

process CELL_BARCODE_CTSS {

    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/ctss/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(bam_sense)
    tuple val(meta), path(bam_antisense)
    tuple val(meta), path(barcodes)      // Filtered cell barcodes

    output:
    tuple val(meta), path("${meta.id}.sense.ctss.bed.gz"),     emit: ctss_sense
    tuple val(meta), path("${meta.id}.antisense.ctss.bed.gz"), emit: ctss_antisense
    tuple val(meta), path("${meta.id}.ctss.bed.gz"),           emit: ctss_bed
    tuple val(meta), path("${meta.id}.ctss.log"),              emit: log

    script:
    """
    #!/bin/bash
    set -euo pipefail

    echo "[CTSS] Generating per-cell CTSS BED for: ${meta.id}" > ${meta.id}.ctss.log

    # ── Sense strand CTSS (5'-most position of forward reads) ──
    # For sense reads: the 5'-most position IS the TSS
    # Extract: chr, TSS_start, TSS_end, cell_barcode, mapq, strand
    samtools view -@ ${task.cpus} -q ${params.min_mapq} ${bam_sense} | \\
        awk '
        {
            # Extract CB (cell barcode) tag
            cb = "";
            for (i=12; i<=NF; i++) {
                if (\$i ~ /^CB:Z:/) { cb = substr(\$i, 6); break }
            }
            if (cb == "") next;
            # 5-prime position for sense reads = leftmost position
            tss_start = \$4 - 1;   # 0-based
            tss_end   = \$4;        # 1-based end
            print \$3 "\\t" tss_start "\\t" tss_end "\\t" cb "\\t" \$5 "\\t+"
        }' | \\
        sort -k4,4 | \\
        # Filter to valid cell barcodes only
        join -1 4 -2 1 - <(sort ${barcodes}) | \\
        awk '{print \$2"\\t"\$3"\\t"\$4"\\t"\$1"\\t"\$5"\\t"\$6}' | \\
        sort -k1,1 -k2,2n | \\
        gzip > ${meta.id}.sense.ctss.bed.gz

    # ── Antisense strand CTSS (5'-most position of reverse reads) ──
    # For antisense reads: the 5' position is the RIGHTMOST alignment pos
    # (because they're on the minus strand)
    samtools view -@ ${task.cpus} -q ${params.min_mapq} ${bam_antisense} | \\
        awk '
        {
            cb = "";
            for (i=12; i<=NF; i++) {
                if (\$i ~ /^CB:Z:/) { cb = substr(\$i, 6); break }
            }
            if (cb == "") next;
            # 5-prime position for antisense = rightmost position
            # Calculate from CIGAR string alignment length
            aln_len = 0;
            cigar = \$6;
            while (match(cigar, /[0-9]+[MDN]/)) {
                aln_len += substr(cigar, RSTART, RLENGTH-1) + 0;
                cigar = substr(cigar, RSTART+RLENGTH);
            }
            tss_end   = \$4 + aln_len;
            tss_start = tss_end - 1;
            print \$3 "\\t" tss_start "\\t" tss_end "\\t" cb "\\t" \$5 "\\t-"
        }' | \\
        sort -k4,4 | \\
        join -1 4 -2 1 - <(sort ${barcodes}) | \\
        awk '{print \$2"\\t"\$3"\\t"\$4"\\t"\$1"\\t"\$5"\\t"\$6}' | \\
        sort -k1,1 -k2,2n | \\
        gzip > ${meta.id}.antisense.ctss.bed.gz

    # ── Merge sense + antisense into single CTSS BED ──────────────
    zcat ${meta.id}.sense.ctss.bed.gz ${meta.id}.antisense.ctss.bed.gz | \\
        sort -k1,1 -k2,2n | \\
        gzip > ${meta.id}.ctss.bed.gz

    SENSE_N=\$(zcat ${meta.id}.sense.ctss.bed.gz | wc -l)
    ANTI_N=\$(zcat ${meta.id}.antisense.ctss.bed.gz | wc -l)
    echo "[CTSS] Sense CTSS records: \$SENSE_N" >> ${meta.id}.ctss.log
    echo "[CTSS] Antisense CTSS records: \$ANTI_N" >> ${meta.id}.ctss.log
    echo "[CTSS] Total CTSS records: \$((\$SENSE_N + \$ANTI_N))" >> ${meta.id}.ctss.log
    """
}
