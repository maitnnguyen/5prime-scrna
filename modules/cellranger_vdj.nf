// ============================================================
//  modules/cellranger_vdj.nf
//  Cell Ranger VDJ — BCR repertoire assembly
//  10x Chromium 5' BCR library
//
//  Notes:
//  - Loaded via HPC module (not Singularity container)
//  - Handles multi-lane FASTQs automatically (L001+L002)
//  - --sample=${meta.id}_BCR matches files named SAMPLE1_BCR_S*_R*
//    while --id=${meta.id} keeps output directory name consistent
//    with the rest of the pipeline (no _BCR suffix in output paths)
//  - Output: clonotypes, contigs, filtered annotations
//  - Can be integrated with Seurat via scRepertoire
//  - Skip with: --skip_cellranger
// ============================================================

process CELLRANGER_VDJ {

    tag "$meta.id (${params.genome})"
    label 'process_high'

    publishDir "${params.outdir}/${params.genome}/cellranger/vdj/${meta.id}", mode: 'copy',
        saveAs: { filename ->
            if (filename =~ /filtered_contig/)    return filename
            if (filename =~ /clonotypes/)         return filename
            if (filename =~ /consensus/)          return filename
            if (filename =~ /web_summary/)        return filename
            if (filename =~ /metrics_summary/)    return filename
            if (filename =~ /airr_rearrangement/) return filename
            if (filename =~ /cloupe/)             return null
            return null
        }

    input:
    // type: 'dir' tells Nextflow this is a directory — required for correct staging
    tuple val(meta), path(fastq_dir, stageAs: 'fastq_dir', type: 'dir')
    path  vdj_ref

    output:
    tuple val(meta), path("${meta.id}/outs/filtered_contig_annotations.csv"), emit: csv
    tuple val(meta), path("${meta.id}/outs/filtered_contig.fasta"),           emit: fasta
    tuple val(meta), path("${meta.id}/outs/filtered_contig.fastq"),           emit: fastq
    tuple val(meta), path("${meta.id}/outs/airr_rearrangement.tsv"),          emit: airr
    tuple val(meta), path("${meta.id}/outs/clonotypes.csv"),                  emit: clonotypes
    tuple val(meta), path("${meta.id}/outs/web_summary.html"),                emit: summary
    tuple val(meta), path("${meta.id}/outs/metrics_summary.csv"),             emit: metrics
    tuple val(meta), path("${meta.id}/outs/consensus_annotations.csv"),       emit: consensus

    script:
    // ────────────────────────────────────────────────────────
    // --id=${meta.id}          : Output dir named by sample ID (no _BCR suffix)
    // --sample=${meta.id}_BCR  : Matches BCR FASTQ prefix e.g. SAMPLE1_BCR_S1_L001_R1_001.fastq.gz
    //                            GEX files (SAMPLE1_*) are automatically excluded
    // --chain=IG               : B-cell receptor only (IGH + IGK + IGL)
    // ────────────────────────────────────────────────────────
    """
    cellranger vdj \\
        --id=${meta.id} \\
        --fastqs=${fastq_dir} \\
        --sample=${meta.id}_BCR \\
        --reference=${vdj_ref} \\
        --chain=IG \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()} \\
        --disable-ui

    echo "[CellRanger VDJ] Done for ${meta.id}"
    """
}
