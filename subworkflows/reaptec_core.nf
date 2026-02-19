// ============================================================
//  subworkflows/reaptec_core.nf
//  Core ReapTEC steps: STAR → UMI dedup → SoftclipG → CTSS
//  Based on: github.com/MurakawaLab/ReapTEC
// ============================================================

include { STAR_ALIGN          } from '../modules/star_align'
include { UMITOOLS_DEDUP      } from '../modules/umitools_dedup'
include { SOFTCLIP_G_FILTER   } from '../modules/softclip_g'
include { CELL_BARCODE_CTSS   } from '../modules/cell_barcode_ctss'
include { CTSS_COUNTS_BIGWIG  } from '../modules/ctss_counts_bigwig'

workflow REAPTEC_CORE {

    take:
    ch_reads      // [ meta, [R1, R2] ] trimmed reads
    ch_star_index // STAR genome index
    ch_gtf        // Gene annotation GTF
    ch_whitelist  // 10x barcode whitelist

    main:

    // ── STEP 1: STARsolo alignment (5' chemistry, strand-aware) ──
    // Key ReapTEC note: Read1 must be >26bp (e.g. 150bp paired-end)
    // to capture the 5' cap signature. STARsolo parameters are
    // specifically tuned for 5' chemistry and UMI/CB extraction.
    STAR_ALIGN(
        ch_reads,
        ch_star_index,
        ch_gtf,
        ch_whitelist
    )

    // ── STEP 2: UMI-tools deduplication ──────────────────────────
    // Remove PCR duplicates while preserving the 5' TSS position.
    // Critical: dedup must be done BEFORE extracting 5' positions
    // to avoid counting the same molecule multiple times.
    UMITOOLS_DEDUP(
        STAR_ALIGN.out.bam,
        STAR_ALIGN.out.bai
    )

    // ── STEP 3: SoftclipG filter (ReapTEC cap signature) ─────────
    // This is the KEY ReapTEC innovation:
    // Extracts reads starting with an unencoded G at the 5' end.
    // The extra G is added by the reverse transcriptase at the
    // 5' cap of the mRNA (template switching) — it marks the
    // EXACT transcription start site (TSS).
    // Sense reads: R1 starts with softclipped G
    // Antisense reads: R1 ends with softclipped C (complementary)
    // This separates genuine 5' cap reads from internal priming.
    SOFTCLIP_G_FILTER(
        UMITOOLS_DEDUP.out.bam,
        UMITOOLS_DEDUP.out.bai
    )

    // ── STEP 4: Cell barcode CTSS BED generation ─────────────────
    // Identifies 5'-end TSS positions (CTSS = CAGE-like TSS)
    // from the softclipG-filtered, STARsolo-aligned,
    // UMI-deduplicated 5' scRNA-seq data.
    // Output: per-cell BED file with columns:
    //   chr  start  end  cell_barcode  count  strand
    CELL_BARCODE_CTSS(
        SOFTCLIP_G_FILTER.out.bam_sense,
        SOFTCLIP_G_FILTER.out.bam_antisense,
        STAR_ALIGN.out.barcodes
    )

    // ── STEP 5: CTSS counts + BigWig normalization ────────────────
    // Counts 5' ends at known promoters and FANTOM5 enhancers.
    // Generates CPM-normalized BigWig files for IGV/UCSC
    // visualization of the TSS landscape.
    CTSS_COUNTS_BIGWIG(
        CELL_BARCODE_CTSS.out.ctss_bed,
        STAR_ALIGN.out.barcodes
    )

    emit:
    bam          = UMITOOLS_DEDUP.out.bam       // Deduplicated BAM
    ctss_bed     = CELL_BARCODE_CTSS.out.ctss_bed  // Per-cell CTSS BED
    bigwig_sense = CTSS_COUNTS_BIGWIG.out.bw_sense  // Sense strand BigWig
    bigwig_anti  = CTSS_COUNTS_BIGWIG.out.bw_anti   // Antisense BigWig
    gene_matrix  = STAR_ALIGN.out.gene_matrix    // STARsolo gene x cell matrix
    barcodes     = STAR_ALIGN.out.barcodes        // Filtered barcodes
    logs         = STAR_ALIGN.out.log
                   .mix(UMITOOLS_DEDUP.out.log)
                   .mix(SOFTCLIP_G_FILTER.out.log)
}
