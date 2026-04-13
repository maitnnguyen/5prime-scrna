// ============================================================
//  subworkflows/reaptec_core.nf
//
//  Step order (matches official ReapTEC log, Murakawa lab 2023):
//    STAR_ALIGN
//      → SOFTCLIP_G_FILTER   (R1 unique mappers with cap-G/C signature)
//        → UMITOOLS_DEDUP    (dedup with CR/UR tags, umi_tools container)
//          → FILTER_BARCODES (filter to valid CB whitelist, HPC samtools)
//            → CTSS_BED      (1-bp TSS map with cell barcodes)
//              ├→ CTSS_COUNTS_BIGWIG   (CPM BigWig + FANTOM QC)
//              └→ BIDIR_ENHANCER_CALL  (bidirectional enhancer calling)
// ============================================================

include { STAR_ALIGN          } from '../modules/star_align'
include { SOFTCLIP_G_FILTER   } from '../modules/softclip_g'
include { UMITOOLS_DEDUP      } from '../modules/umitools_dedup'
include { CELL_BARCODE_FILTER     } from '../modules/cell_barcode_filter'
include { CTSS_BED            } from '../modules/ctss_bed'
include { CTSS_COUNTS_BIGWIG  } from '../modules/ctss_counts_bigwig'
include { BIDIR_ENHANCER_CALL } from '../modules/bidir_enhancer_call'

workflow REAPTEC_CORE {

    take:
    ch_star_input   // [ meta, r1_files, r2_files, whitelist ]
    ch_star_index   // STAR genome index (value channel)
    ch_gtf          // Gene annotation GTF (value channel)
    ch_chrom_sizes  // genome.chrom.sizes (value channel)
    ch_fantom_promo // FANTOM5 promoter BED — [] for T2T
    ch_fantom_enhan // FANTOM-NET enhancer BED — [] for T2T
    ch_mask_bed     // Masking BED for enhancer calling (value channel)

    main:

    // Unpack input — split reads and whitelist into separate channels
    ch_reads     = ch_star_input.map { meta, r1, r2, wl -> [ meta, r1, r2 ] }
    ch_whitelist = ch_star_input.map { meta, r1, r2, wl -> [ meta, wl ]     }

    // ── STEP 1: STARsolo alignment ────────────────────────────
    STAR_ALIGN(
        ch_reads,
        ch_star_index,
        ch_gtf,
        ch_whitelist
    )

    // ── STEP 2: SoftclipG filter ──────────────────────────────
    // Extract R1 unique mappers with cap-G/C signature
    // Must happen BEFORE dedup — only TSS reads are deduplicated
    SOFTCLIP_G_FILTER(
        STAR_ALIGN.out.bam,
        STAR_ALIGN.out.bai
    )

    // ── STEP 3: UMI deduplication ────────────────────────────
    // CR/UR tags per official ReapTEC protocol
    // umi_tools container — no samtools inside
    UMITOOLS_DEDUP(
        SOFTCLIP_G_FILTER.out.bam
            .join( SOFTCLIP_G_FILTER.out.bai, by: 0 )
    )

    // ── STEP 4: Filter to valid cell barcodes ────────────────
    // HPC samtools -D CB for exact tag matching against whitelist
    //ch_filter_input = UMITOOLS_DEDUP.out.bam
    //    .join( ch_whitelist, by: 0 )

    CELL_BARCODE_FILTER(UMITOOLS_DEDUP.out.bam, ch_whitelist)

    // ── STEP 5: CTSS BED generation ──────────────────────────
    ch_ctss_input = CELL_BARCODE_FILTER.out.bam
        .join( CELL_BARCODE_FILTER.out.bai, by: 0 )

    CTSS_BED(ch_ctss_input)

    // ── STEP 6: CPM BigWig + FANTOM quantification ───────────
    CTSS_COUNTS_BIGWIG(
        CTSS_BED.out.bed,
        ch_chrom_sizes,
        ch_fantom_promo,
        ch_fantom_enhan
    )

    // ── STEP 7: Bidirectional enhancer calling ───────────────
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
