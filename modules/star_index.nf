// ============================================================
//  modules/star_index.nf
//  Build STAR genome index using GENCODE annotation
//  Reference: STARsolo_STARindex_GENCODE41_PRI_221224_logFile
// ============================================================

process STAR_INDEX {

    label 'process_high'
    publishDir "${params.outdir}/star_index", mode: 'copy'

    input:
    path fasta    // Genome FASTA (e.g. GRCh38.primary_assembly.genome.fa.gz)
    path gtf      // Gene annotation GTF (e.g. GENCODE v41 comprehensive)

    output:
    path "star_index/", emit: index
    path "star_index.log", emit: log

    script:
    """
    mkdir -p star_index

    echo "[STAR Index] Building genome index..." > star_index.log
    echo "[STAR Index] Genome: ${fasta}" >> star_index.log
    echo "[STAR Index] GTF: ${gtf}" >> star_index.log

    # Decompress if gzipped
    FASTA_FILE=${fasta}
    if [[ ${fasta} == *.gz ]]; then
        gunzip -c ${fasta} > genome.fa
        FASTA_FILE=genome.fa
    fi

    GTF_FILE=${gtf}
    if [[ ${gtf} == *.gz ]]; then
        gunzip -c ${gtf} > annotation.gtf
        GTF_FILE=annotation.gtf
    fi

    STAR \\
        --runMode genomeGenerate \\
        --runThreadN ${task.cpus} \\
        --genomeDir star_index/ \\
        --genomeFastaFiles \$FASTA_FILE \\
        --sjdbGTFfile \$GTF_FILE \\
        --sjdbOverhang 149 \\
        2>> star_index.log

    echo "[STAR Index] Index built successfully." >> star_index.log
    """
}
