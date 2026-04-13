// ============================================================
//  subworkflows/reaptec_core.nf
//
//  Step order (matches official ReapTEC log, Murakawa lab 2023):
//    STAR_ALIGN
//      → SOFTCLIP_G_FILTER   (R1 unique mappers with cap-G/C signature)
//        → UMITOOLS_DEDUP    (dedup within cap-G reads using CB/UB tags)
//          → CTSS_BED        (1-bp TSS map with cell barcodes)
//            ├→ CTSS_COUNTS_BIGWIG   (CPM BigWig + FANTOM QC)
//            └→ BIDIR_ENHANCER_CALL  (bidirectional enhancer calling)
// ============================================================

include { STAR_ALIGN          } from '../modules/star_align'
include { SOFTCLIP_G_FILTER   } from '../modules/softclip_g'
include { UMITOOLS_DEDUP      } from '../modules/umitools_dedup'
include { CELL_BARCODE_FILTER } from '../modules/cell_barcode_filter'
include { CTSS_BED            } from '../modules/ctss_bed'
include { CTSS_COUNTS_BIGWIG  } from '../modules/ctss_counts_bigwig'
include { BIDIR_ENHANCER_CALL } from '../modules/bidir_enhancer_call'

workflow REAPTEC_CORE {

    take:
    ch_star_input   // [ meta, [r1_files], [r2_files], whitelist ]
                    // Built in main.nf — r1/r2 are lists to support multi-lane
    ch_star_index   // STAR genome index (value channel)
    ch_gtf          // Gene annotation GTF (value channel)
    ch_chrom_sizes  // genome.chrom.sizes (value channel)
    ch_fantom_promo // FANTOM5 promoter BED — [] for T2T (optional)
    ch_fantom_enhan // FANTOM-NET enhancer BED — [] for T2T (optional)
    ch_mask_bed     // Masking BED for enhancer calling (value channel)

    main:

    // ── Unpack input channel ──────────────────────────────────
    // ch_star_input carries [ meta, r1_files, r2_files, whitelist ]
    // We need reads and whitelist as separate channels for STAR_ALIGN
    ch_reads     = ch_star_input.map { meta, r1, r2, wl -> [ meta, r1, r2 ] }
    ch_whitelist = ch_star_input.map { meta, r1, r2, wl -> [ meta, wl ]     }

    // ── STEP 1: STARsolo alignment ────────────────────────────
    STAR_ALIGN(
        ch_reads,       // [ meta, [r1_files], [r2_files] ]
        ch_star_index,
        ch_gtf,
        ch_whitelist    // [ meta, whitelist ] — joined by meta in STAR_ALIGN
    )

    // ── STEP 2: SoftclipG filter ──────────────────────────────
    // Extracts R1-only uniquely mapped reads with cap-G/C signature.
    // Must happen BEFORE dedup — only cap-signature reads are deduplicated.
    SOFTCLIP_G_FILTER(
        STAR_ALIGN.out.bam,
        STAR_ALIGN.out.bai
    )

    // ── STEP 3: UMI deduplication ─────────────────────────────
    // Join bam + bai + whitelist by meta.id before calling.
    // Guarantees correct whitelist is paired with correct sample BAM.
    ch_dedup_input = SOFTCLIP_G_FILTER.out.bam
        .join( SOFTCLIP_G_FILTER.out.bai, by: 0 )
        .join( ch_whitelist,              by: 0 )
        // → [ meta, bam, bai, whitelist ]

    UMITOOLS_DEDUP(ch_dedup_input)

    // ── STEP 3.5: Filter to valid cell barcodes ────────────────
    ch_filter_input = UMITOOLS_DEDUP.out.bam
        .join( ch_whitelist, by: 0 )

    CELL_BARCODE_FILTER(ch_filter_input)

    // ── STEP 4: CTSS BED generation ───────────────────────────
    // Join bam + bai into single tuple so meta is never lost
    ch_ctss_input = CELL_BARCODE_FILTER.out.bam
        .join( CELL_BARCODE_FILTER.out.bai, by: 0 )
        // → [ meta, bam, bai ]

    CTSS_BED(ch_ctss_input)

    // ── STEP 5: CPM BigWig + FANTOM quantification ────────────
    // fantom_promo and fantom_enhan are [] for T2T —
    // CTSS_COUNTS_BIGWIG skips quantification when they are empty
    CTSS_COUNTS_BIGWIG(
        CTSS_BED.out.bed,
        ch_chrom_sizes,
        ch_fantom_promo,
        ch_fantom_enhan
    )

    // ── STEP 6: Bidirectional enhancer calling ────────────────
    BIDIR_ENHANCER_CALL(
        CTSS_BED.out.bed,
        ch_mask_bed
    )

    emit:
    ctss_bed    = CTSS_BED.out.bed
    bigwig_fwd  = CTSS_COUNTS_BIGWIG.out.bw_fwd
    bigwig_rev  = CTSS_COUNTS_BIGWIG.out.bw_rev
    enhancers   = BIDIR_ENHANCER_CALL.out.bed
    tpm_matrix  = BIDIR_ENHANCER_CALL.out.tpm_matrix
    prom_counts = CTSS_COUNTS_BIGWIG.out.prom_counts
    gene_matrix = STAR_ALIGN.out.gene_matrix
}
