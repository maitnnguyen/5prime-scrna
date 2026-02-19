// ============================================================
//  subworkflows/scafe_workflow.nf
//  SCAFE subworkflow:
//    Per-sample solo → aggregate across samples →
//    cross-validate with ReapTEC btcEnh BED
// ============================================================

include { SCAFE_SOLO            } from '../modules/scafe'
include { SCAFE_AGGREGATE       } from '../modules/scafe'
include { SCAFE_COUNT_AGGREGATE } from '../modules/scafe'

workflow SCAFE_WORKFLOW {

    take:
    ch_bam        // [ meta, bam ] — deduplicated BAM from STARsolo
    ch_bai        // [ meta, bai ]
    ch_barcodes   // [ meta, barcodes.tsv.gz ]

    main:

    // ── STEP 1: Run SCAFE solo per sample ────────────────────
    // Runs the full SCAFE pipeline per sample:
    //   bam_to_ctss → remove_strand_invader → cluster
    //   → filter (logistic regression) → annotate → count
    SCAFE_SOLO(
        ch_bam,
        ch_bai,
        ch_barcodes,
        params.scafe_genome
    )

    // ── STEP 2: Aggregate CTSS across ALL samples ─────────────
    // Pool CTSS BEDs from all samples to define a COMMON set of
    // tCREs — essential for multi-sample CVID vs HC comparison
    // so tCRE IDs are portable between samples
    ch_all_ctss = SCAFE_SOLO.out.ctss_bed
        .map { meta, ctss -> ctss }
        .collect()

    SCAFE_AGGREGATE(
        ch_all_ctss,
        params.scafe_genome
    )

    // ── STEP 3: Count UMI per sample against aggregated tCREs ─
    // Re-count per sample using the aggregated (common) tCRE set
    SCAFE_COUNT_AGGREGATE(
        SCAFE_SOLO.out.ctss_bed,
        SCAFE_AGGREGATE.out.tCRE_bed,
        params.scafe_genome
    )

    emit:
    tCRE_bed       = SCAFE_AGGREGATE.out.tCRE_bed        // Aggregated tCRE BED (consensus)
    tss_clusters   = SCAFE_SOLO.out.tss_clusters          // Per-sample filtered TSS clusters
    count_matrices = SCAFE_COUNT_AGGREGATE.out.matrix      // Per-sample UMI matrices
    ctss_beds      = SCAFE_SOLO.out.ctss_bed               // Per-sample CTSS BEDs
    logs           = SCAFE_SOLO.out.log
                     .mix(SCAFE_AGGREGATE.out.log)
}
