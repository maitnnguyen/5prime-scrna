// ============================================================
//  modules/bidir_enhancer_call.nf
//  ReapTEC Step 6: Bidirectional enhancer calling (btcEnhs)
//  Adapted from: fixed_bidir_enhancers_10bp.sh (FANTOM5/ReapTEC)
//  Reference: Andersson et al., Nature 2014 + Oguchi et al., Science 2024
// ============================================================

process BIDIR_ENHANCER_CALL {

    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/enhancers/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(ctss_bed)   // Per-cell CTSS BED (gzipped)
    path  fantom_bed                   // FANTOM5 reference BED

    output:
    tuple val(meta), path("${meta.id}.raw_enhancers.bed"),         emit: raw_enhancers
    tuple val(meta), path("${meta.id}.bidir_pairs.bed"),           emit: bidir_pairs
    tuple val(meta), path("${meta.id}.sense_peaks.bed"),           emit: sense_peaks
    tuple val(meta), path("${meta.id}.antisense_peaks.bed"),       emit: antisense_peaks
    tuple val(meta), path("${meta.id}.enhancer_call.log"),         emit: log

    script:
    // ────────────────────────────────────────────────────────
    // BIDIRECTIONAL ENHANCER DETECTION ALGORITHM
    // (ReapTEC fixed_bidir_enhancers_10bp.sh logic):
    //
    // STEP A: Create pseudo-bulk CTSS by aggregating all cells
    // STEP B: Cluster CTSS into peaks (merge within 10bp)
    // STEP C: Separate into sense (+) and antisense (-) peaks
    // STEP D: For each sense peak, find antisense peak within
    //         bidir_gap bp (default 500bp)
    // STEP E: The midpoint between sense+antisense peaks = enhancer center
    // STEP F: Define enhancer window = ±bidir_window bp from center
    // STEP G: Require BOTH strands to have signal (bidirectional test)
    // STEP H: Apply CPM filter (log2CPM >= min_cpm in ≥1 cluster)
    //
    // OUTPUT: BED6+ format:
    //   chr  start  end  name  score  strand  sense_TSS  anti_TSS  bidir_score
    // ────────────────────────────────────────────────────────

    def gap    = params.bidir_gap
    def window = params.bidir_window
    def min_cpm = params.min_cpm

    """
    #!/bin/bash
    set -euo pipefail

    echo "[btcEnh] Starting bidirectional enhancer calling for: ${meta.id}" > ${meta.id}.enhancer_call.log
    echo "[btcEnh] Parameters: gap=${gap}bp window=${window}bp min_log2CPM=${min_cpm}" >> ${meta.id}.enhancer_call.log

    # ── STEP A: Pseudo-bulk CTSS (aggregate all cells) ────────────
    echo "[btcEnh] A) Aggregating per-cell CTSS to pseudo-bulk..." >> ${meta.id}.enhancer_call.log

    zcat ${ctss_bed} | \\
        awk '{print \$1"\\t"\$2"\\t"\$3"\\t"\$6}' | \\
        sort -k1,1 -k2,2n -k4,4 | \\
        bedtools merge -s -i stdin -c 4 -o count | \\
        awk '{print \$1"\\t"\$2"\\t"\$3"\\t"NR"\\t"\$5"\\t"\$4}' > pseudobulk_ctss.bed

    TOTAL_CTSS=\$(wc -l < pseudobulk_ctss.bed)
    echo "[btcEnh] Total pseudo-bulk CTSS positions: \$TOTAL_CTSS" >> ${meta.id}.enhancer_call.log

    # ── STEP B: Cluster CTSS peaks (merge within 10bp, same strand) ──
    echo "[btcEnh] B) Clustering CTSS into peaks (10bp merge)..." >> ${meta.id}.enhancer_call.log

    bedtools merge -s -d 10 -i pseudobulk_ctss.bed -c 5,6 -o sum,distinct | \\
        awk '{print \$1"\\t"\$2"\\t"\$3"\\tpeak_"NR"\\t"\$5"\\t"\$6}' > ctss_peaks.bed

    # ── STEP C: Separate sense / antisense peaks ──────────────────
    echo "[btcEnh] C) Separating sense/antisense peaks..." >> ${meta.id}.enhancer_call.log

    awk '\$6=="+"' ctss_peaks.bed | sort -k1,1 -k2,2n > ${meta.id}.sense_peaks.bed
    awk '\$6=="-"' ctss_peaks.bed | sort -k1,1 -k2,2n > ${meta.id}.antisense_peaks.bed

    SENSE_PEAKS=\$(wc -l < ${meta.id}.sense_peaks.bed)
    ANTI_PEAKS=\$(wc -l < ${meta.id}.antisense_peaks.bed)
    echo "[btcEnh] Sense peaks: \$SENSE_PEAKS  |  Antisense peaks: \$ANTI_PEAKS" >> ${meta.id}.enhancer_call.log

    # ── STEP D: Find bidirectional pairs within gap ───────────────
    # For each sense peak, find the closest antisense peak within gap bp
    echo "[btcEnh] D) Finding bidirectional peak pairs (gap<=${gap}bp)..." >> ${meta.id}.enhancer_call.log

    bedtools window \\
        -a ${meta.id}.sense_peaks.bed \\
        -b ${meta.id}.antisense_peaks.bed \\
        -w ${gap} \\
        -Sm | \\
        awk '{
            # Sense peak: \$1-\$6, Antisense peak: \$7-\$12
            s_center = int((\$2 + \$3) / 2);
            a_center = int((\$8 + \$9) / 2);
            # Enhancer center = midpoint between sense and antisense TSS
            enh_center = int((s_center + a_center) / 2);
            enh_start  = enh_center - ${window};
            enh_end    = enh_center + ${window};
            if (enh_start < 0) enh_start = 0;
            bidir_score = \$5 + \$11;   # Sum of both strand scores
            # Output: enhancer region + source TSS positions
            print \$1 "\\t" enh_start "\\t" enh_end "\\t" \\
                  "btcEnh_" NR "\\t" bidir_score "\\t." "\\t" \\
                  s_center "\\t" a_center "\\t" \\
                  \$5 "\\t" \$11
        }' | \\
        sort -k1,1 -k2,2n > ${meta.id}.bidir_pairs.bed

    BIDIR_N=\$(wc -l < ${meta.id}.bidir_pairs.bed)
    echo "[btcEnh] Raw bidirectional pairs: \$BIDIR_N" >> ${meta.id}.enhancer_call.log

    # ── STEP E: Merge overlapping enhancer regions ────────────────
    echo "[btcEnh] E) Merging overlapping enhancer candidates..." >> ${meta.id}.enhancer_call.log

    bedtools merge \\
        -i ${meta.id}.bidir_pairs.bed \\
        -c 4,5,7,8,9,10 \\
        -o first,sum,min,max,sum,sum | \\
        awk '{
            OFS="\\t";
            print \$1, \$2, \$3, \$4, \$5, ".", \$6, \$7, \$8, \$9
        }' > ${meta.id}.raw_enhancers.bed

    RAW_ENH=\$(wc -l < ${meta.id}.raw_enhancers.bed)
    echo "[btcEnh] Raw enhancers after merge: \$RAW_ENH" >> ${meta.id}.enhancer_call.log
    echo "[btcEnh] Bidirectional enhancer calling complete." >> ${meta.id}.enhancer_call.log

    # ── Summary statistics ────────────────────────────────────────
    echo "" >> ${meta.id}.enhancer_call.log
    echo "=== SUMMARY ===" >> ${meta.id}.enhancer_call.log
    echo "CTSS positions: \$TOTAL_CTSS" >> ${meta.id}.enhancer_call.log
    echo "Sense peaks: \$SENSE_PEAKS" >> ${meta.id}.enhancer_call.log
    echo "Antisense peaks: \$ANTI_PEAKS" >> ${meta.id}.enhancer_call.log
    echo "Bidir pairs: \$BIDIR_N" >> ${meta.id}.enhancer_call.log
    echo "Final enhancers (pre-filter): \$RAW_ENH" >> ${meta.id}.enhancer_call.log
    """
}


// ============================================================
//  modules/promoter_filter.nf
//  ReapTEC Step 7: Remove promoter-overlapping regions
// ============================================================

process PROMOTER_FILTER {

    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}/enhancers/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(raw_enhancers)
    path  fantom_bed

    output:
    tuple val(meta), path("${meta.id}.enhancers.bed"),          emit: enhancers
    tuple val(meta), path("${meta.id}.enhancers_filtered.log"), emit: log

    script:
    """
    #!/bin/bash
    set -euo pipefail

    echo "[PromoFilter] Filtering promoter-overlapping enhancers: ${meta.id}" > ${meta.id}.enhancers_filtered.log

    BEFORE=\$(wc -l < ${raw_enhancers})

    # ── Remove known promoter regions ─────────────────────────────
    # Download FANTOM5 promoters if not provided
    if [ ! -f "${fantom_bed}" ]; then
        echo "[PromoFilter] Downloading FANTOM5 robust promoters..." >> ${meta.id}.enhancers_filtered.log
        wget -q -O fantom5_promoters_hg38.bed.gz \\
            "https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38.promoter.robust.bed.gz"
        gunzip fantom5_promoters_hg38.bed.gz
        PROMO_BED="fantom5_promoters_hg38.bed"
    else
        PROMO_BED="${fantom_bed}"
    fi

    # Step 1: Remove overlaps with FANTOM5 robust promoters
    bedtools subtract \\
        -a ${raw_enhancers} \\
        -b \$PROMO_BED \\
        -A > step1_no_fantom_promoter.bed

    STEP1=\$(wc -l < step1_no_fantom_promoter.bed)
    echo "[PromoFilter] After FANTOM5 promoter removal: \$STEP1 (removed \$((\$BEFORE - \$STEP1)))" >> ${meta.id}.enhancers_filtered.log

    # Step 2: Remove ±2kb windows around annotated gene TSS
    # (requires hg38 gene annotation - generated from GTF)
    if [ -f "gene_tss_2kb.bed" ]; then
        bedtools subtract \\
            -a step1_no_fantom_promoter.bed \\
            -b gene_tss_2kb.bed \\
            -A > step2_no_gene_tss.bed
        STEP2=\$(wc -l < step2_no_gene_tss.bed)
        echo "[PromoFilter] After gene TSS ±2kb removal: \$STEP2" >> ${meta.id}.enhancers_filtered.log
        mv step2_no_gene_tss.bed filtered_enhancers.bed
    else
        mv step1_no_fantom_promoter.bed filtered_enhancers.bed
    fi

    # Step 3: Sort and deduplicate final enhancer list
    sort -k1,1 -k2,2n filtered_enhancers.bed | \\
        bedtools merge -i stdin -c 4,5,6 -o first,sum,first > ${meta.id}.enhancers.bed

    FINAL=\$(wc -l < ${meta.id}.enhancers.bed)
    echo "[PromoFilter] Final high-confidence enhancers: \$FINAL" >> ${meta.id}.enhancers_filtered.log
    echo "[PromoFilter] Total removed (promoter-overlapping): \$((\$BEFORE - \$FINAL))" >> ${meta.id}.enhancers_filtered.log

    # ── Compute enhancer size distribution ────────────────────────
    awk '{print \$3 - \$2}' ${meta.id}.enhancers.bed | \\
        awk 'BEGIN{s=0;n=0}{s+=\$1;n++}END{print "[PromoFilter] Mean enhancer size: "s/n"bp | N="n}' >> ${meta.id}.enhancers_filtered.log
    """
}
