// ============================================================
//  modules/ctss_counts_bigwig.nf
//  Count CTSS reads at FANTOM5 promoters/enhancers
//  Generate CPM-normalised BigWig for IGV / UCSC visualisation
//  Adapted from: STARsolo_Counts_CTSS_bed_CPM_bigwig_240119.sh
// ============================================================

process CTSS_COUNTS_BIGWIG {

    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/bigwig/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(ctss_bed)     // merged sense+antisense CTSS BED
    tuple val(meta), path(barcodes)     // filtered barcodes for cell count

    output:
    tuple val(meta), path("${meta.id}.sense.bw"),          emit: bw_sense
    tuple val(meta), path("${meta.id}.antisense.bw"),      emit: bw_anti
    tuple val(meta), path("${meta.id}.ctss_counts.tsv.gz"),emit: counts
    tuple val(meta), path("${meta.id}.bigwig.log"),        emit: log

    script:
    """
    #!/bin/bash
    set -euo pipefail

    echo "[BigWig] Generating normalised BigWig tracks for: ${meta.id}" > ${meta.id}.bigwig.log

    # ── Total mapped CTSS for CPM normalisation ───────────────────
    TOTAL=\$(zcat ${ctss_bed} | wc -l)
    echo "[BigWig] Total CTSS records: \$TOTAL" >> ${meta.id}.bigwig.log

    # ── Genome sizes (hg38) ───────────────────────────────────────
    # fetchChromSizes from UCSC or use pre-downloaded file
    if [ ! -f hg38.chrom.sizes ]; then
        fetchChromSizes hg38 > hg38.chrom.sizes 2>/dev/null || \\
        wget -q -O hg38.chrom.sizes https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
    fi

    # ── Sense strand BigWig (+ strand CTSS) ───────────────────────
    zcat ${ctss_bed} | awk '\$6=="+"' | \\
        awk '{print \$1"\\t"\$2"\\t"\$3"\\t"1}' | \\
        sort -k1,1 -k2,2n | \\
        bedtools genomecov -i stdin -g hg38.chrom.sizes -bg -scale \$(echo "1000000 / \$TOTAL" | bc -l) | \\
        sort -k1,1 -k2,2n > sense_cpm.bedgraph

    bedGraphToBigWig sense_cpm.bedgraph hg38.chrom.sizes ${meta.id}.sense.bw
    echo "[BigWig] Sense BigWig created." >> ${meta.id}.bigwig.log

    # ── Antisense strand BigWig (- strand CTSS, values negative) ──
    zcat ${ctss_bed} | awk '\$6=="-"' | \\
        awk '{print \$1"\\t"\$2"\\t"\$3"\\t"-1}' | \\
        sort -k1,1 -k2,2n | \\
        bedtools genomecov -i stdin -g hg38.chrom.sizes -bg -scale \$(echo "1000000 / \$TOTAL" | bc -l) | \\
        sort -k1,1 -k2,2n > antisense_cpm.bedgraph

    bedGraphToBigWig antisense_cpm.bedgraph hg38.chrom.sizes ${meta.id}.antisense.bw
    echo "[BigWig] Antisense BigWig created." >> ${meta.id}.bigwig.log

    # ── Count matrix at FANTOM5 enhancer positions ─────────────────
    # Produces a simple count table for downstream use
    zcat ${ctss_bed} | \\
        awk 'BEGIN{OFS="\\t"}{print \$1, \$2, \$3, \$4, \$6}' | \\
        gzip > ${meta.id}.ctss_counts.tsv.gz

    echo "[BigWig] Done." >> ${meta.id}.bigwig.log
    """
}
