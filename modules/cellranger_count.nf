// ============================================================
//  modules/cellranger_count.nf
//  Cell Ranger count — GEX gene expression quantification
//  10x Chromium 5' scRNA-seq
//
//  Notes:
//  - Loaded via HPC module (not Singularity container)
//  - Handles multi-lane FASTQs automatically (L001+L002)
//  - --sample must match the FASTQ filename prefix exactly
//    e.g. SAMPLE1_S1_L001_R1_001.fastq.gz → --sample=SAMPLE1
//  - Output BAM used by: SCAFE, CamoTSS
//  - Output matrix used by: Seurat clustering, EXTRACT_WHITELIST
//  - Skip with: --skip_cellranger
// ============================================================

process CELLRANGER_COUNT {

    tag "$meta.id (${params.genome})"
    label 'process_high'

    publishDir "${params.outdir}/${params.genome}/cellranger/${meta.id}", mode: 'copy',
        saveAs: { filename ->
            if (filename =~ /possorted_genome_bam/)  return "bam/$filename"
            if (filename =~ /filtered_feature_bc/)   return "matrix/$filename"
            if (filename =~ /raw_feature_bc/)        return null
            if (filename =~ /web_summary/)           return filename
            if (filename =~ /metrics_summary/)       return filename
            if (filename =~ /molecule_info/)         return filename
            if (filename =~ /cloupe/)                return null
            return null
        }

    input:
    // type: 'dir' tells Nextflow this is a directory — required for correct staging
    tuple val(meta), path(fastq_dir, stageAs: 'fastq_dir', type: 'dir')
    path  cellranger_ref

    output:
    tuple val(meta), path("${meta.id}/outs/possorted_genome_bam.bam"),                   emit: bam
    tuple val(meta), path("${meta.id}/outs/possorted_genome_bam.bam.bai"),               emit: bai
    tuple val(meta), path("${meta.id}/outs/filtered_feature_bc_matrix/"),                emit: matrix
    tuple val(meta), path("${meta.id}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"), emit: barcodes
    tuple val(meta), path("${meta.id}/outs/web_summary.html"),                           emit: summary
    tuple val(meta), path("${meta.id}/outs/metrics_summary.csv"),                        emit: metrics
    tuple val(meta), path("${meta.id}/outs/molecule_info.h5"),                           emit: molecule_info

    script:
    // ────────────────────────────────────────────────────────
    // --chemistry SC5P-PE   : 5' paired-end chemistry
    //                         DO NOT use 'fiveprime' — strips
    //                         TSO sequence needed for SCAFE/CamoTSS
    // --include-introns     : Include intronic reads
    // --create-bam true     : Required for SCAFE and CamoTSS
    // --sample=${meta.id}   : Matches FASTQ prefix e.g. SAMPLE1_S1_L001_R1_001.fastq.gz
    //                         BCR files (SAMPLE1_BCR_*) are automatically excluded
    //                         because they don't match this prefix
    // ────────────────────────────────────────────────────────
    """
    cellranger count \\
        --id=${meta.id} \\
        --fastqs=${fastq_dir} \\
        --sample=${meta.id} \\
        --transcriptome=${cellranger_ref} \\
        --chemistry=SC5P-PE \\
        --include-introns=true \\
        --create-bam=true \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()} \\
        --disable-ui

    echo "[CellRanger count] Done for ${meta.id}"
    """
}
