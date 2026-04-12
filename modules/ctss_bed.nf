// ============================================================
//  modules/ctss_bed.nf
//  Generate CTSS BED file with Cell Barcodes for 10x 5' scRNA-seq
//  Adapted from: STARsolo_Cell_barcode_CTSS_bed_20221123.sh
//
//  Input BAM comes from UMITOOLS_DEDUP (after softclip + dedup).
//  Only cap-G reads with valid barcodes reach this step.
// ============================================================

process CTSS_BED {

    tag "$meta.id"
    label 'process_medium'

    publishDir "${params.outdir}/${params.genome}/reaptec/ctss", mode: 'copy'

    input:
    // bam and bai arrive as a single joined tuple — meta is never lost
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.CTSS.bed"), emit: bed

    script:
    """
    # Extract 5' TSS position + Cell Barcode directly from BAM
    # CB tag = corrected barcode (written by STARsolo after whitelist match)
    samtools view ${bam} | awk 'BEGIN{OFS="\\t"} {
        barcode = "unknown"
        for (i = 12; i <= NF; i++) {
            if (\$i ~ /^CB:Z:/) { barcode = substr(\$i, 6) }
        }

        # Reverse strand (-): TSS is at the END of the alignment
        if (and(\$2, 16)) {
            print \$3, \$4 + length(\$10) - 1, \$4 + length(\$10), barcode, ".", "-"
        } else {
        # Forward strand (+): TSS is at the START of the alignment
            print \$3, \$4 - 1, \$4, barcode, ".", "+"
        }
    }' | sort -k1,1 -k2,2n -k6,6 > temp_tss.bed

    # Group by coordinate + barcode + strand, count reads per position per cell
    bedtools groupby -i temp_tss.bed -g 1,2,3,4,6 -c 1 -o count | \\
        awk 'BEGIN{OFS="\\t"} {
            # Output: Chr | Start | End | Barcode | Count | Strand
            print \$1, \$2, \$3, \$4, \$6, \$5
        }' > ${meta.id}.CTSS.bed

    rm temp_tss.bed
    """
}
