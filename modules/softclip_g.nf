// ============================================================
//  modules/softclip_g.nf
//  ReapTEC Step 3: Cap signature (softclipG) filter
//  Adapted from: MurakawaLab/ReapTEC/bin/STARsolo_SoftclipG_221224.sh
// ============================================================

process SOFTCLIP_G_FILTER {

    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/softclip_g/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(bam)
    tuple val(meta), path(bai)

    output:
    tuple val(meta), path("SoftclipG_${meta.id}.bam"), emit: bam
    path "SoftclipG_${meta.id}.bam.bai"              , emit: bai
    tuple val(meta), path("${meta.id}.softclipG.log"),     emit: log

    script:
    // ────────────────────────────────────────────────────────
    // THE REAPTEC CAP SIGNATURE CONCEPT:
    //
    // When reverse transcriptase reaches the 5' cap of an mRNA,
    // it template-switches and adds non-templated C residues.
    // These become G residues in the cDNA/read.
    //
    // In 5' scRNA-seq (10x Chromium 5' kit):
    //   → Sense reads: The 5' G appears as a SOFTCLIPPED G at
    //     the START of the aligned read (CIGAR: NNS...M)
    //     where N = softclipped bases at 5' end
    //
    //   → Antisense reads: The complementary C appears as a
    //     SOFTCLIPPED C at the END of the aligned read
    //     (CIGAR: M...NNS at 3' end)
    //
    // By extracting ONLY reads with this cap signature, ReapTEC:
    //   1. Pinpoints exact TSS positions at nucleotide resolution
    //   2. Separates genuine 5' cap reads from internal priming
    //   3. Enables bidirectional enhancer detection
    //
    // Position of cap G in read:
    //   Next GEM: position 40 of R1 (after 28bp CB+UMI + 11bp clip)
    //   GEM-X:    position 42 of R1 (change BASE and CIGAR pattern)
    // ────────────────────────────────────────────────────────

    // clip_val is the base position (e.g., 40)
    // cigar_val is the string match (e.g., 40S)
    def clip_val  = params.clip5p + 1 
    def cigar_val = params.clip5p + 1 

    """
    # Initialize log
    echo "[ReapTEC SoftclipG] Sample: ${meta.id}" > ${meta.id}.softclipG.log
    echo "[ReapTEC SoftclipG] Target Position: ${clip_val}" >> ${meta.id}.softclipG.log

    # 1. Extract Unique Read 1s 
    # -f 64: Read 1 only
    # -q 255: Unique mappers only
    samtools view -@ ${task.cpus} -h -f 64 -q 255 ${bam} > unique_r1.sam

    # 2. Filter for Softclip G (Forward) and C (Reverse)
    # We use a subshell ( ... ) to stream both filters into one sort command
    (
        samtools view -H unique_r1.sam; # Keep the BAM header
        
        # Sense/Forward: starts with softclip, base at position is G
        samtools view -F 16 unique_r1.sam | awk -F '\\t' 'BEGIN {OFS="\\t"} {
            BASE = substr(\$10, ${clip_val}, 1);
            if (\$6 ~ /^${cigar_val}S/ && (BASE == "G" || BASE == "g")) {print \$0}
        }';
        
        # Antisense/Reverse: ends with softclip, base at calculated position is C
        samtools view -f 16 unique_r1.sam | awk -F '\\t' 'BEGIN {OFS="\\t"} {
            ALT = substr(\$10, length(\$10)-(${clip_val}-1), 1);
            if (\$6 ~ /${cigar_val}S\$/ && (ALT == "C" || ALT == "c")) {print \$0}
        }'
    ) | samtools sort -@ ${task.cpus} -o SoftclipG_${meta.id}.bam

    # 3. Index the final ReapTEC-ready BAM
    samtools index SoftclipG_${meta.id}.bam
    
    # Cleanup
    rm unique_r1.sam
    echo "[ReapTEC SoftclipG] Filter complete." >> ${meta.id}.softclipG.log
    """
}
