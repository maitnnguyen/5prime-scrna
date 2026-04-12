// ============================================================
//  modules/star_align.nf
//  STARsolo alignment — ReapTEC 5' scRNA-seq parameters
//
//  Key parameters validated against official ReapTEC log
//  (Murakawa lab, April 2023):
//    --soloStrand Reverse     (5' chemistry — NOT Forward)
//    --soloUMIlen 10          (NOT 12; GEM-X: use 12)
//    --clip5pNbases 39 0      (16 CB + 10 UMI + 13 linker; GEM-X: 41)
//    --soloBarcodeMate 1      (barcode on R1)
//    Read order: R1 R2        (NOT R2 R1)
//    --soloStrand Reverse     (confirmed in official log — NOT Forward)
// ============================================================

process STAR_ALIGN {

    tag "$meta.id (${params.genome})"
    label 'process_high'

    publishDir "${params.outdir}/${params.genome}/star/${meta.id}", mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith('.log'))      "logs/$filename"
            else if (filename.endsWith('.bam')) "bam/$filename"
            else if (filename =~ /Solo\.out/)  "solo_out/$filename"
            else null
        }

    input:
    tuple val(meta), path(r1), path(r2)  // R1: CB+UMI+TSS; R2: cDNA
                                          // r1/r2 can be lists for multi-lane
    path  star_index                      // STAR genome index directory
    path  gtf                             // Gene annotation GTF
    tuple val(meta), path(whitelist)      // Per-sample whitelist from EXTRACT_WHITELIST

    output:
    tuple val(meta), path("${meta.id}.Aligned.sortedByCoord.out.bam"),     emit: bam
    tuple val(meta), path("${meta.id}.Aligned.sortedByCoord.out.bam.bai"), emit: bai
    tuple val(meta), path("${meta.id}.Solo.out/"),                          emit: solo_out
    tuple val(meta), path("${meta.id}.Solo.out/Gene/filtered/"),            emit: gene_matrix
    tuple val(meta), path("${meta.id}.Log.final.out"),                      emit: log
    tuple val(meta), path("${meta.id}.SJ.out.tab"),                         emit: junctions

    script:
    // For multi-lane: STAR accepts comma-separated file lists
    def r1_files = (r1 instanceof List ? r1 : [r1]).join(',')
    def r2_files = (r2 instanceof List ? r2 : [r2]).join(',')

    """
    # Allow STAR to open many files (required for BAM sorting on HPC)
    ulimit -n 10000

    STAR \\
        --runMode alignReads \\
        --runThreadN ${task.cpus} \\
        --genomeDir ${star_index} \\
        --sjdbGTFfile ${gtf} \\
        --readFilesIn ${r1_files} ${r2_files} \\
        --readFilesCommand zcat \\
        --soloType CB_UMI_Simple \\
        --soloCBwhitelist ${whitelist} \\
        --soloBarcodeMate 1 \\
        --soloCBstart  ${params.cb_start} \\
        --soloCBlen    ${params.cb_len} \\
        --soloUMIstart ${params.umi_start} \\
        --soloUMIlen   ${params.umi_len} \\
        --soloStrand   ${params.solo_strand} \\
        --clip5pNbases ${params.clip5p} 0 \\
        --soloFeatures Gene GeneFull \\
        --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \\
        --soloUMIdedup 1MM_Directional_UMItools \\
        --alignEndsType Local \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \\
        --outSAMstrandField intronMotif \\
        --outFilterIntronMotifs RemoveNoncanonicalUnannotated \\
        --outFilterMultimapNmax 1 \\
        --outFilterMismatchNoverLmax 0.05 \\
        --alignSJDBoverhangMin 1 \\
        --alignIntronMax 1000000 \\
        --outFileNamePrefix ${meta.id}. \\
        --limitBAMsortRAM ${task.memory.toBytes() - 2000000000}

    samtools index -@ ${task.cpus} ${meta.id}.Aligned.sortedByCoord.out.bam
    """
}
