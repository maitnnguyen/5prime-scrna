// ============================================================
//  modules/star_align.nf
//  STARsolo alignment — ReapTEC 5' scRNA-seq parameters
// ============================================================

process STAR_ALIGN {

    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}/star/${meta.id}", mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith('.log'))       "logs/$filename"
            else if (filename.endsWith('.bam'))  "bam/$filename"
            else if (filename =~ /Solo.out/)     "solo_out/$filename"
            else null
        }

    input:
    tuple val(meta), path(reads)    // [meta, [R1_trimmed, R2_trimmed]]
    path  star_index                // STAR genome index directory
    path  gtf                       // Gene annotation GTF
    path  whitelist                 // 10x barcode whitelist

    output:
    tuple val(meta), path("${meta.id}.Aligned.sortedByCoord.out.bam"), emit: bam
    tuple val(meta), path("${meta.id}.Aligned.sortedByCoord.out.bam.bai"), emit: bai
    tuple val(meta), path("${meta.id}.Solo.out/"), emit: solo_out
    tuple val(meta), path("${meta.id}.Solo.out/Gene/filtered/"), emit: gene_matrix
    tuple val(meta), path("${meta.id}.Solo.out/Gene/filtered/barcodes.tsv"), emit: barcodes
    tuple val(meta), path("${meta.id}.Log.final.out"), emit: log
    tuple val(meta), path("${meta.id}.SJ.out.tab"), emit: junctions

    script:
    // ────────────────────────────────────────────────────────
    // ReapTEC critical parameters (from STARsolo_ReapTEC_logfile):
    //
    // --soloType CB_UMI_Simple      : 10x Chromium barcode+UMI
    // --soloCBwhitelist              : 10x v3 3M barcode list
    // --soloCBstart/len              : CB = bases 1-16 of R1
    // --soloUMIstart/len             : UMI = bases 17-28 of R1
    // --soloStrand Forward           : 5' chemistry is FORWARD
    // --clip5pNbases 39 0            : Remove 39bp from R1 5' end
    //                                  (28bp CB+UMI + 11bp linker)
    //                                  This exposes the cDNA for
    //                                  TSS identification
    //                                  GEM-X: change 39 → 41
    // --alignEndsType Local          : Allows softclipping at 5' end
    //                                  Essential for detecting the
    //                                  unencoded G cap signature!
    // --outSAMattributes CB UB ...   : Keep cell/UMI tags in BAM
    // ────────────────────────────────────────────────────────
    def r1 = reads[0]   // Read 1: CB + UMI (28bp) + TSS
    def r2 = reads[1]   // Read 2: cDNA (150bp)

    """
    # ── STARsolo 5' scRNA-seq alignment (ReapTEC parameters) ──
    STAR \\
        --runMode alignReads \\
        --runThreadN ${task.cpus} \\
        --genomeDir ${star_index} \\
        --sjdbGTFfile ${gtf} \\
        --readFilesIn ${r2} ${r1} \\
        --readFilesCommand zcat \\
        \\
        --soloType CB_UMI_Simple \\
        --soloCBwhitelist ${whitelist} \\
        --soloCBstart ${params.cb_start} \\
        --soloCBlen ${params.cb_len} \\
        --soloUMIstart ${params.umi_start} \\
        --soloUMIlen ${params.umi_len} \\
        --soloStrand ${params.solo_strand} \\
        --soloFeatures Gene Velocyto \\
        --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \\
        --soloUMIfiltering MultiGeneUMI_CR \\
        --soloCellFilter EmptyDrops_CR \\
        \\
        --clip5pNbases ${params.clip5p} 0 \\
        --alignEndsType Local \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMattributes NH HI AS NM MD CB UB \\
        --outSAMstrandField intronMotif \\
        --outFilterIntronMotifs RemoveNoncanonicalUnannotated \\
        --outFilterMultimapNmax 1 \\
        --outFilterMismatchNoverLmax 0.05 \\
        --alignSJDBoverhangMin 1 \\
        --alignIntronMax 1000000 \\
        \\
        --outFileNamePrefix ${meta.id}. \\
        --outSAMheaderHD @HD VN:1.4 SO:coordinate \\
        --limitBAMsortRAM ${task.memory.toBytes() - 2000000000}

    # ── Index BAM ──────────────────────────────────────────────
    samtools index -@ ${task.cpus} ${meta.id}.Aligned.sortedByCoord.out.bam

    # ── Add CB tag to reads for downstream filtering ───────────
    # Note: STARsolo already adds CB/UB tags; verify presence
    samtools view -H ${meta.id}.Aligned.sortedByCoord.out.bam | \\
        grep -q "CB:Z" || echo "WARNING: CB tag not found in BAM"
    """
}
