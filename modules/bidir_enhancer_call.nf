// ============================================================
//  modules/bidir_enhancer_call.nf
//  ReapTEC Step 6: Bidirectional enhancer calling (btcEnhs)
//  Adapted from: fixed_bidir_enhancers_10bp.sh (FANTOM5/ReapTEC)
//  Reference: Andersson et al., Nature 2014 + Oguchi et al., Science 2024
// ============================================================

process BIDIR_ENHANCER_CALL {

    tag "$meta.id"
    label 'process_medium'
    
    // Organizes by genome (T2T vs hg38)
    publishDir "${params.outdir}/${params.genome}/reaptec/enhancers", mode: 'copy'

    input:
    tuple val(meta), path(ctss_bed) // The .CTSS.bed file from previous module
    path mask_bed                  // The Gencode (hg38) or T2T-specific mask BED

    output:
    tuple val(meta), path("${meta.id}_enhancers/"), emit: dir
    path "${meta.id}_enhancers/*.enhancers.bed",    emit: bed
    path "${meta.id}_enhancers/*.expression.tpm.matrix", emit: tpm_matrix
    path "${meta.id}_enhancers/*.directionality.txt",    emit: directionality

    script:
    """
    # 1. Prepare the bedlist file required by the Andersson/Murakawa script
    # The script expects a file containing absolute paths to input BEDs
    echo "\$(readlink -f ${ctss_bed})" > bedlist.txt

    # 2. Create the output directory
    mkdir -p ${meta.id}_enhancers

    # 3. Execute the official ReapTEC enhancer call script
    # This script is located in your project's bin/ folder (automatically in PATH)
    # -f: path to the bedlist
    # -m: the masking file (crucial for T2T vs hg38)
    # -s: prefix/stub for output files
    # -o: output directory
    
    fixed_bidir_enhancers_10bp \\
        -f bedlist.txt \\
        -m ${mask_bed} \\
        -s ${meta.id}_ \\
        -o ${meta.id}_enhancers/
        
    # Note: If the script fails with 'command not found', 
    # check that chmod +x was run on all bin/ files.
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
