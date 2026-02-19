// ============================================================
//  subworkflows/enhancer_analysis.nf
//  Bidirectional enhancer calling + B cell downstream analysis
// ============================================================

include { BIDIR_ENHANCER_CALL   } from '../modules/bidir_enhancer_call'
include { PROMOTER_FILTER       } from '../modules/promoter_filter'
include { SEURAT_CLUSTERING     } from '../modules/seurat_clustering'
include { PSEUDOBULK_ENHANCER   } from '../modules/pseudobulk_enhancer'

workflow ENHANCER_ANALYSIS {

    take:
    ch_ctss       // Per-cell CTSS BED files
    ch_fantom     // FANTOM5 promoter/enhancer reference
    ch_bam        // Deduplicated BAM (for pseudobulk)
    ch_gex        // Gene expression matrix (STARsolo)
    ch_barcodes   // Filtered cell barcodes

    main:

    // ── STEP 6: Bidirectional enhancer calling ────────────────────
    // Identifies btcEnhs (bidirectionally transcribed candidate
    // enhancers) from the CTSS data.
    // Algorithm (adapted from FANTOM5 fixed_bidir_enhancers_10bp.sh):
    //   1. Merge CTSS peaks within 10bp on same strand
    //   2. Find sense/antisense peak pairs within bidir_gap bp
    //   3. Require signal on BOTH strands (bidirectional)
    //   4. Remove pairs overlapping annotated promoters (FANTOM5)
    //   5. Apply log2CPM >= min_cpm filter for robust enhancers
    BIDIR_ENHANCER_CALL(
        ch_ctss,
        ch_fantom
    )

    // ── STEP 7: Promoter masking ──────────────────────────────────
    // Strictly filter enhancer candidates that overlap with:
    //   - Known gene promoters (±2kb from TSS)
    //   - FANTOM5 robust promoters
    //   - Repeat elements (optional)
    // Remaining peaks = high-confidence bidirectional enhancers
    PROMOTER_FILTER(
        BIDIR_ENHANCER_CALL.out.raw_enhancers,
        ch_fantom
    )

    // ── STEP 8: B cell clustering (Seurat) ───────────────────────
    // Cluster cells using gene expression matrix.
    // Annotate B cell subtypes:
    //   - Transitional B cells
    //   - Naïve B cells       (CD27- IgD+)
    //   - Switched Memory     (CD27+ IgD-)  ← depleted in CVID
    //   - Unswitched Memory   (CD27+ IgD+)
    //   - GC B cells          (BCL6+ AICDA+)
    //   - Plasmablasts        (PRDM1+ XBP1+)
    SEURAT_CLUSTERING(
        ch_gex,
        ch_barcodes,
        PROMOTER_FILTER.out.enhancers,
        ch_ctss
    )

    // ── STEP 9: Pseudobulk enhancer analysis ─────────────────────
    // Aggregate CTSS counts per B cell cluster.
    // Run DESeq2 differential enhancer analysis: CVID vs Control
    // per B cell subtype.
    // Outputs per-cluster enhancer activity matrices for
    // downstream TF motif and GRN analysis.
    PSEUDOBULK_ENHANCER(
        PROMOTER_FILTER.out.enhancers,
        SEURAT_CLUSTERING.out.cluster_assignments,
        ch_ctss,
        ch_bam
    )

    emit:
    bidir_enhancers     = PROMOTER_FILTER.out.enhancers
    ctss_matrix         = PSEUDOBULK_ENHANCER.out.count_matrix
    diff_enhancers      = PSEUDOBULK_ENHANCER.out.diff_results
    cluster_assignments = SEURAT_CLUSTERING.out.cluster_assignments
    seurat_obj          = SEURAT_CLUSTERING.out.seurat_rds
    logs                = BIDIR_ENHANCER_CALL.out.log
                          .mix(SEURAT_CLUSTERING.out.log)
                          .mix(PSEUDOBULK_ENHANCER.out.log)
}
