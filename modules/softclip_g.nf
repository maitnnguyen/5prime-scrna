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
    tuple val(meta), path("${meta.id}.sense.bam"),     emit: bam_sense
    tuple val(meta), path("${meta.id}.sense.bam.bai"), emit: bai_sense
    tuple val(meta), path("${meta.id}.antisense.bam"),     emit: bam_antisense
    tuple val(meta), path("${meta.id}.antisense.bam.bai"), emit: bai_antisense
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

    def cap_pos    = params.clip5p + 1           // e.g. 40 for Next GEM
    def cigar_pos  = params.clip5p               // e.g. 39 for softclip CIGAR
    def cigar_alt  = params.clip5p               // antisense end

    """
    #!/bin/bash
    set -euo pipefail

    echo "[ReapTEC SoftclipG] Processing: ${meta.id}" > ${meta.id}.softclipG.log
    echo "[ReapTEC SoftclipG] Cap position in read: ${cap_pos}" >> ${meta.id}.softclipG.log
    echo "[ReapTEC SoftclipG] Expected CIGAR pattern: ${cigar_pos}S[0-9]+M..." >> ${meta.id}.softclipG.log

    # ── Extract SENSE reads (forward strand, cap G at 5' start) ──
    # Filter: read is on forward strand (-F 16) AND
    #         the base at position cap_pos is 'G' AND
    #         the CIGAR starts with the expected softclip pattern
    echo "[ReapTEC] Extracting sense reads with 5' cap G..." >> ${meta.id}.softclipG.log

    samtools view -@ ${task.cpus} -F 16 ${bam} | \\
        awk -v pos=${cap_pos} -v clip=${cigar_pos} '
        {
            BASE = substr(\$10, pos, 1);
            if (\$6 ~ "^" clip "S[0-9]" && BASE == "G") { print \$0 }
        }' | \\
        cat <(samtools view -H ${bam}) - | \\
        samtools sort -@ ${task.cpus} -o ${meta.id}.sense.bam

    samtools index -@ ${task.cpus} ${meta.id}.sense.bam

    SENSE_COUNT=\$(samtools view -c ${meta.id}.sense.bam)
    echo "[ReapTEC] Sense reads with cap G: \$SENSE_COUNT" >> ${meta.id}.softclipG.log

    # ── Extract ANTISENSE reads (reverse strand, cap C at 3' end) ──
    # The complementary C at the 3' end of antisense reads
    echo "[ReapTEC] Extracting antisense reads with 3' cap C..." >> ${meta.id}.softclipG.log

    samtools view -@ ${task.cpus} -f 16 ${bam} | \\
        awk -v clip=${cigar_alt} '
        {
            ALT = substr(\$10, length(\$10) - clip + 1, 1);
            if (\$6 ~ "[0-9]M" clip "S\$" && ALT == "C") { print \$0 }
        }' | \\
        cat <(samtools view -H ${bam}) - | \\
        samtools sort -@ ${task.cpus} -o ${meta.id}.antisense.bam

    samtools index -@ ${task.cpus} ${meta.id}.antisense.bam

    ANTI_COUNT=\$(samtools view -c ${meta.id}.antisense.bam)
    echo "[ReapTEC] Antisense reads with cap C: \$ANTI_COUNT" >> ${meta.id}.softclipG.log

    TOTAL=\$(samtools view -c ${bam})
    echo "[ReapTEC] Total input reads: \$TOTAL" >> ${meta.id}.softclipG.log
    echo "[ReapTEC] Cap signature retention rate: \$(echo "scale=2; (\$SENSE_COUNT + \$ANTI_COUNT) / \$TOTAL * 100" | bc)%" >> ${meta.id}.softclipG.log

    echo "[ReapTEC SoftclipG] Done." >> ${meta.id}.softclipG.log
    """
}
