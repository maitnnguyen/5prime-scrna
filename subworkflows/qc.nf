// ============================================================
//  subworkflows/qc.nf
//  QC subworkflow: FastQC → Trim Galore → FastQC (post-trim)
// ============================================================

include { FASTQC    } from '../modules/fastqc'
include { TRIMGALORE } from '../modules/trimgalore'

workflow QC_SUBWORKFLOW {

    take:
    ch_reads   // [ meta, [R1, R2] ]

    main:

    // Pre-trim QC
    FASTQC(ch_reads)

    // Adapter trimming
    TRIMGALORE(ch_reads)

    // Collect all logs for MultiQC
    ch_logs = FASTQC.out.zip
        .mix(TRIMGALORE.out.log)

    emit:
    reads = TRIMGALORE.out.reads   // Trimmed reads for alignment
    logs  = ch_logs                 // QC logs for MultiQC
}
