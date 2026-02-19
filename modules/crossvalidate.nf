// ============================================================
//  modules/crossvalidate.nf
//  Cross-validate ReapTEC btcEnh BED vs SCAFE tCRE BED
//  Produces tiered consensus BED for scTSS input
//
//  Confidence tiers:
//    TIER 1 — Validated by BOTH ReapTEC (cap G) + SCAFE (logistic regression)
//    TIER 2 — ReapTEC only (bidirectional, cap G validated)
//    TIER 3 — SCAFE only (ML validated, not bidirectional)
// ============================================================

process CROSSVALIDATE_ENHANCERS {

    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/crossvalidate/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(reaptec_bed)  // ReapTEC btcEnh BED (from BIDIR_ENHANCER_CALL)
    tuple val(meta), path(scafe_bed)    // SCAFE aggregated tCRE BED (distal only)

    output:
    tuple val(meta), path("${meta.id}.tier1_both.bed"),       emit: tier1       // Both tools
    tuple val(meta), path("${meta.id}.tier2_reaptec.bed"),    emit: tier2       // ReapTEC only
    tuple val(meta), path("${meta.id}.tier3_scafe.bed"),      emit: tier3       // SCAFE only
    tuple val(meta), path("${meta.id}.consensus.bed"),         emit: consensus   // All tiers merged
    tuple val(meta), path("${meta.id}.crossval.log"),          emit: log

    script:
    """
    #!/bin/bash
    set -euo pipefail

    echo "[CrossVal] Cross-validating ReapTEC btcEnh vs SCAFE tCRE for: ${meta.id}" > ${meta.id}.crossval.log

    # ── Extract SCAFE distal tCREs (enhancers only) ───────────────
    # SCAFE annotates tCREs as 'proximal' (promoter) or 'distal' (enhancer)
    # We only want distal tCREs for comparison with ReapTEC btcEnhs
    if [[ "${scafe_bed}" == *.gz ]]; then
        zcat ${scafe_bed} | awk '\$6=="distal"' > scafe_distal.bed
    else
        awk '\$6=="distal"' ${scafe_bed} > scafe_distal.bed
    fi

    SCAFE_DISTAL=\$(wc -l < scafe_distal.bed)
    REAPTEC_N=\$(wc -l < ${reaptec_bed})
    echo "[CrossVal] ReapTEC btcEnhs: \$REAPTEC_N" >> ${meta.id}.crossval.log
    echo "[CrossVal] SCAFE distal tCREs: \$SCAFE_DISTAL" >> ${meta.id}.crossval.log

    # ── TIER 1: Overlapping in BOTH tools ─────────────────────────
    # These are highest-confidence enhancers:
    #   - Biologically validated (cap G signature from ReapTEC)
    #   - Statistically validated (logistic regression from SCAFE)
    bedtools intersect \\
        -a ${reaptec_bed} \\
        -b scafe_distal.bed \\
        -u \\
        | awk 'BEGIN{OFS="\\t"}{print \$0, "TIER1_BOTH"}' \\
        | sort -k1,1 -k2,2n > ${meta.id}.tier1_both.bed

    TIER1=\$(wc -l < ${meta.id}.tier1_both.bed)
    echo "[CrossVal] TIER 1 (validated by BOTH): \$TIER1" >> ${meta.id}.crossval.log

    # ── TIER 2: ReapTEC only (not in SCAFE) ───────────────────────
    # Bidirectional cap-signature enhancers not captured by SCAFE
    # (may be low-expression or missed by SCAFE clustering)
    bedtools intersect \\
        -a ${reaptec_bed} \\
        -b scafe_distal.bed \\
        -v \\
        | awk 'BEGIN{OFS="\\t"}{print \$0, "TIER2_REAPTEC"}' \\
        | sort -k1,1 -k2,2n > ${meta.id}.tier2_reaptec.bed

    TIER2=\$(wc -l < ${meta.id}.tier2_reaptec.bed)
    echo "[CrossVal] TIER 2 (ReapTEC only):       \$TIER2" >> ${meta.id}.crossval.log

    # ── TIER 3: SCAFE only (not in ReapTEC) ───────────────────────
    # ML-validated tCREs not confirmed as bidirectional
    # (may be unidirectional active enhancers or unannotated promoters)
    bedtools intersect \\
        -a scafe_distal.bed \\
        -b ${reaptec_bed} \\
        -v \\
        | awk 'BEGIN{OFS="\\t"}{print \$0, "TIER3_SCAFE"}' \\
        | sort -k1,1 -k2,2n > ${meta.id}.tier3_scafe.bed

    TIER3=\$(wc -l < ${meta.id}.tier3_scafe.bed)
    echo "[CrossVal] TIER 3 (SCAFE only):         \$TIER3" >> ${meta.id}.crossval.log

    # ── CONSENSUS: Merge all tiers ────────────────────────────────
    # Unified BED for downstream scTSS / Seurat input
    # Tier annotation in column 7 (or last column)
    cat ${meta.id}.tier1_both.bed \\
        ${meta.id}.tier2_reaptec.bed \\
        ${meta.id}.tier3_scafe.bed \\
        | sort -k1,1 -k2,2n \\
        | bedtools merge -i stdin -c 4,5,7 -o first,sum,first \\
        > ${meta.id}.consensus.bed

    CONSENSUS=\$(wc -l < ${meta.id}.consensus.bed)
    echo "[CrossVal] Consensus BED (all tiers): \$CONSENSUS" >> ${meta.id}.crossval.log

    # ── Summary stats ─────────────────────────────────────────────
    echo "" >> ${meta.id}.crossval.log
    echo "=== CROSS-VALIDATION SUMMARY ===" >> ${meta.id}.crossval.log
    echo "ReapTEC btcEnhs:           \$REAPTEC_N" >> ${meta.id}.crossval.log
    echo "SCAFE distal tCREs:        \$SCAFE_DISTAL" >> ${meta.id}.crossval.log
    echo "TIER 1 (BOTH):             \$TIER1  ← highest confidence" >> ${meta.id}.crossval.log
    echo "TIER 2 (ReapTEC only):     \$TIER2" >> ${meta.id}.crossval.log
    echo "TIER 3 (SCAFE only):       \$TIER3" >> ${meta.id}.crossval.log
    echo "Consensus total:           \$CONSENSUS" >> ${meta.id}.crossval.log
    PCT_OVERLAP=\$(echo "scale=1; \$TIER1 / \$REAPTEC_N * 100" | bc 2>/dev/null || echo "N/A")
    echo "ReapTEC validated by SCAFE: \$PCT_OVERLAP%" >> ${meta.id}.crossval.log
    echo "[CrossVal] Done." >> ${meta.id}.crossval.log
    """
}
