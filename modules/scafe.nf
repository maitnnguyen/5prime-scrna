// ============================================================
//  modules/scafe.nf
//  SCAFE — Single Cell Analysis of Five-prime Ends
//  github.com/chung-lab/SCAFE
//  Moody & Kouno et al., Bioinformatics 2022
//
//  SCAFE workflow (solo mode):
//    BAM → bam_to_ctss → remove_strand_invader → cluster
//        → filter (logistic regression, AUC>0.98)
//        → annotate → count → tCRE BED + UMI matrix
// ============================================================

process SCAFE_SOLO {

    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}/scafe/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(bam)          // Sorted, CB/UB-tagged BAM (from STARsolo or Cell Ranger)
    tuple val(meta), path(bai)          // BAM index
    tuple val(meta), path(barcodes)     // Filtered cell barcodes (barcodes.tsv.gz)
    val   genome                        // SCAFE genome name, e.g. hg38.gencode_v41

    output:
    tuple val(meta), path("${meta.id}/count/*/matrix/"),            emit: tCRE_matrix
    tuple val(meta), path("${meta.id}/annotate/*/bed/*.tCRE.bed.gz"), emit: tCRE_bed
    tuple val(meta), path("${meta.id}/filter/*/bed/*.tssCluster.default.filtered.bed.gz"), emit: tss_clusters
    tuple val(meta), path("${meta.id}/bam_to_ctss/*/bed/*.collapse.ctss.bed.gz"), emit: ctss_bed
    tuple val(meta), path("${meta.id}.scafe.log"),                  emit: log

    script:
    // ────────────────────────────────────────────────────────
    // SCAFE workflow.sc.solo steps (from chung-lab/SCAFE):
    //
    // 1. bam_to_ctss     — extract 5' capped TSS positions from BAM
    //                      auto-detects TS oligo mode (match/trim/skip)
    //                      generates *.ctss.bed.gz and *.unencoded_G.ctss.bed.gz
    //
    // 2. remove_strand_invader — align TS oligo seq to upstream region
    //                            remove strand invasion artifacts
    //
    // 3. cluster         — Paraclu clustering of CTSS into TSS clusters
    //                      min_cluster_count=5, min_summit_count=3
    //
    // 4. filter          — Multiple logistic regression classifier (AUC>0.98)
    //                      uses pre-trained model or --training_signal_path ATAC
    //                      filters genuine TSS from exon-painting artifacts
    //
    // 5. annotate        — merge nearby TSS clusters into tCREs
    //                      classify as proximal (promoter) or distal (enhancer)
    //                      proximity_slop_rng=500bp from annotated gene TSS
    //
    // 6. count           — count UMI per tCRE per cell barcode
    //                      output: tCRE × cell UMI matrix (Market Exchange format)
    //
    // KEY SCAFE NOTE:
    //   BAM must have CB:Z and UB:Z tags (from Cell Ranger OR STARsolo)
    //   TS oligo must be intact on Read 1 (do NOT use --chemistry fiveprime
    //   in Cell Ranger, use SC5P-PE instead)
    //   If using STARsolo BAM: verify CB/UB tags are present before running
    // ────────────────────────────────────────────────────────

    def atac_flag  = params.scafe_atac_bw ? "--training_signal_path ${params.scafe_atac_bw} --testing_signal_path ${params.scafe_atac_bw}" : ""
    def logr_model = params.scafe_glm_model ? "--usr_glm_model_path ${params.scafe_glm_model}" : ""

    """
    #!/bin/bash
    set -euo pipefail

    echo "[SCAFE] Starting solo workflow for: ${meta.id}" > ${meta.id}.scafe.log
    echo "[SCAFE] BAM: ${bam}" >> ${meta.id}.scafe.log
    echo "[SCAFE] Genome: ${genome}" >> ${meta.id}.scafe.log
    echo "[SCAFE] Threads: ${task.cpus}" >> ${meta.id}.scafe.log

    # ── Verify CB/UB tags present in BAM ─────────────────────────
    # SCAFE requires CB:Z and UB:Z tags from Cell Ranger or STARsolo
    CB_CHECK=\$(samtools view ${bam} | head -1000 | grep -c "CB:Z:" || true)
    if [ "\$CB_CHECK" -eq 0 ]; then
        echo "[SCAFE] WARNING: CB:Z tag not detected in BAM - SCAFE requires CB/UB tags" >> ${meta.id}.scafe.log
        echo "[SCAFE] Ensure BAM was produced with STARsolo (--outSAMattributes CB UB) or Cell Ranger SC5P-PE" >> ${meta.id}.scafe.log
    else
        echo "[SCAFE] CB:Z tags confirmed present in BAM" >> ${meta.id}.scafe.log
    fi

    # ── Decompress barcodes if needed ────────────────────────────
    if [[ "${barcodes}" == *.gz ]]; then
        gunzip -c ${barcodes} > barcodes_plain.tsv
        BARCODES_PATH=barcodes_plain.tsv
    else
        BARCODES_PATH=${barcodes}
    fi

    # ── Run SCAFE solo workflow ───────────────────────────────────
    scafe.workflow.sc.solo \\
        --overwrite=yes \\
        --run_bam_path=${bam} \\
        --run_cellbarcode_path=\$BARCODES_PATH \\
        --genome=${genome} \\
        --run_tag=${meta.id} \\
        --run_outDir=${meta.id}/ \\
        --max_thread=${task.cpus} \\
        --detect_TS_oligo=auto \\
        ${atac_flag} \\
        ${logr_model} \\
        2>> ${meta.id}.scafe.log

    echo "[SCAFE] Solo workflow complete." >> ${meta.id}.scafe.log

    # ── Log key output file sizes ─────────────────────────────────
    echo "" >> ${meta.id}.scafe.log
    echo "=== KEY OUTPUTS ===" >> ${meta.id}.scafe.log

    # tCRE BED
    TCRE_BED=\$(find ${meta.id}/annotate/ -name "*.tCRE.bed.gz" 2>/dev/null | head -1)
    if [ -n "\$TCRE_BED" ]; then
        N_TCRE=\$(zcat \$TCRE_BED | wc -l)
        echo "[SCAFE] tCREs defined: \$N_TCRE" >> ${meta.id}.scafe.log
        # Count proximal vs distal
        PROXIMAL=\$(zcat \$TCRE_BED | awk '\$6=="proximal"' | wc -l || true)
        DISTAL=\$(zcat \$TCRE_BED | awk '\$6=="distal"' | wc -l || true)
        echo "[SCAFE] Proximal tCREs (promoters): \$PROXIMAL" >> ${meta.id}.scafe.log
        echo "[SCAFE] Distal tCREs (enhancers):   \$DISTAL" >> ${meta.id}.scafe.log
    fi

    # TSS clusters after filter
    FILTERED_TSS=\$(find ${meta.id}/filter/ -name "*.default.filtered.bed.gz" 2>/dev/null | head -1)
    if [ -n "\$FILTERED_TSS" ]; then
        N_TSS=\$(zcat \$FILTERED_TSS | wc -l)
        echo "[SCAFE] Genuine TSS clusters (post logistic filter): \$N_TSS" >> ${meta.id}.scafe.log
    fi

    # UMI count matrix
    MATRIX_DIR=\$(find ${meta.id}/count/ -name "matrix.mtx.gz" -printf "%h\\n" 2>/dev/null | head -1)
    if [ -n "\$MATRIX_DIR" ]; then
        echo "[SCAFE] UMI count matrix: \$MATRIX_DIR" >> ${meta.id}.scafe.log
    fi

    echo "[SCAFE] Done for: ${meta.id}" >> ${meta.id}.scafe.log
    """
}


// ============================================================
//  process SCAFE_AGGREGATE
//  Aggregate CTSS across multiple samples → unified tCRE set
//  Recommended for multi-sample studies (CVID vs HC)
//  Input: *.ctss.bed.gz files from SCAFE_SOLO bam_to_ctss step
// ============================================================

process SCAFE_AGGREGATE {

    label 'process_high'
    publishDir "${params.outdir}/scafe/aggregate", mode: 'copy'

    input:
    path ctss_beds    // All *.collapse.ctss.bed.gz files collected from all samples
    val  genome

    output:
    path "aggregate/annotate/*/bed/*.tCRE.bed.gz",                          emit: tCRE_bed
    path "aggregate/filter/*/bed/*.tssCluster.default.filtered.bed.gz",     emit: tss_clusters
    path "aggregate.scafe.log",                                              emit: log

    script:
    // ────────────────────────────────────────────────────────
    // Aggregation is recommended for multi-sample studies:
    //   1. Pool CTSS from all samples
    //   2. Re-run cluster → filter → annotate on the pooled data
    //   3. Produces a COMMON set of tCREs for all samples
    //      → portable tCRE IDs across CVID patients and HC
    //
    // This is the recommended SCAFE approach for differential
    // analysis between conditions (Moody & Kouno et al., 2022)
    // ────────────────────────────────────────────────────────

    """
    #!/bin/bash
    set -euo pipefail

    echo "[SCAFE aggregate] Aggregating \$(ls *.ctss.bed.gz | wc -l) CTSS files..." > aggregate.scafe.log

    # List all ctss bed files into an input file
    ls *.collapse.ctss.bed.gz | tr '\\n' ',' | sed 's/,\$//' > ctss_list.txt
    echo "[SCAFE aggregate] CTSS files: \$(cat ctss_list.txt)" >> aggregate.scafe.log

    scafe.workflow.cm.aggregate \\
        --overwrite=yes \\
        --ctss_bed_paths=\$(cat ctss_list.txt) \\
        --genome=${genome} \\
        --run_tag=aggregate \\
        --run_outDir=aggregate/ \\
        --max_thread=${task.cpus} \\
        2>> aggregate.scafe.log

    N_TCRE=\$(find aggregate/annotate/ -name "*.tCRE.bed.gz" | xargs zcat 2>/dev/null | wc -l || echo "0")
    echo "[SCAFE aggregate] Total tCREs in aggregated set: \$N_TCRE" >> aggregate.scafe.log
    echo "[SCAFE aggregate] Done." >> aggregate.scafe.log
    """
}


// ============================================================
//  process SCAFE_COUNT_AGGREGATE
//  Count UMI per aggregated tCRE per cell for each sample
//  Run after SCAFE_AGGREGATE to get per-sample count matrices
//  using the common aggregated tCRE definitions
// ============================================================

process SCAFE_COUNT_AGGREGATE {

    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/scafe/count_aggregate/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(ctss_bed)   // Per-sample ctss.bed.gz from SCAFE_SOLO
    path  tCRE_bed                    // Aggregated tCRE BED from SCAFE_AGGREGATE
    val   genome

    output:
    tuple val(meta), path("${meta.id}/count/*/matrix/"), emit: matrix
    tuple val(meta), path("${meta.id}.count.log"),       emit: log

    script:
    """
    #!/bin/bash
    set -euo pipefail

    echo "[SCAFE count] Counting UMI for: ${meta.id} against aggregated tCREs" > ${meta.id}.count.log

    scafe.tool.sc.count \\
        --overwrite=yes \\
        --ctss_bed_path=${ctss_bed} \\
        --tCRE_bed_path=${tCRE_bed} \\
        --genome=${genome} \\
        --outputPrefix=${meta.id} \\
        --outDir=${meta.id}/ \\
        2>> ${meta.id}.count.log

    echo "[SCAFE count] Done for ${meta.id}" >> ${meta.id}.count.log
    """
}
