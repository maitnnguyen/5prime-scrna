#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * ============================================================
 *  Cell Ranger Dual Pipeline: GEX + VDJ (BCR)
 *  - cellranger count for gene expression
 *  - cellranger vdj for B-cell receptors
 * ============================================================
 */

// Parameters
params.fastq_dir      = null                // Root dir containing all FASTQs
params.transcriptome  = /home/arkku/group/ics/tools/refdata-gex-GRCh38-2024-A  // For cellranger count (GEX)
params.vdj_reference  = /home/arkku/group/ics/tools/cellranger/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0  // For cellranger vdj (BCR)
params.chemistry      = "auto"
params.outdir         = "${projectDir}/results"

// Container paths (update with your actual filenames)
params.container_cache = '/home/arkku/group/ics/tools/singularity_cache'
params.fastqc_container     = "${params.container_cache}/fastqc.sif"        // UPDATE THIS
params.multiqc_container    = "${params.container_cache}/multiqc.sif"       // UPDATE THIS
params.cellranger_container = "${params.container_cache}/cellranger.sif"

// FastQC
process FASTQC {
    tag "${sample_id}"
    label 'process_low'
    publishDir "${params.outdir}/fastqc/${sample_id}", mode: 'copy'

    container params.fastqc_container

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.html"), emit: html
    tuple val(sample_id), path("*.zip"),  emit: zip
    path  "*.zip",                         emit: reports

    script:
    """
    fastqc \\
        --threads ${task.cpus} \\
        --outdir . \\
        ${reads}
    """
}

// MultiQC
process MULTIQC {
    label 'process_low'
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    container params.multiqc_container

    input:
    path reports

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data/",       emit: data

    script:
    """
    multiqc . --filename multiqc_report.html --force
    """
}

// Cell Ranger count (Gene Expression)
process CELLRANGER_COUNT {
    tag "${sample_id}"
    label 'process_high'
    publishDir "${params.outdir}/cellranger_count/${sample_id}", mode: 'copy'

    container params.cellranger_container

    input:
    tuple val(sample_id), path(fastq_dir)
    path  transcriptome

    output:
    tuple val(sample_id), path("${sample_id}/outs/filtered_feature_bc_matrix/"),
        emit: filtered_matrix
    tuple val(sample_id), path("${sample_id}/outs/web_summary.html"),
        emit: web_summary
    tuple val(sample_id), path("${sample_id}/outs/metrics_summary.csv"),
        emit: metrics

    script:
    def chem_arg = params.chemistry != 'auto' ? "--chemistry=${params.chemistry}" : ""
    """
    cellranger count \\
        --id=${sample_id} \\
        --transcriptome=${transcriptome} \\
        --fastqs=${fastq_dir} \\
        --sample=${sample_id} \\
        ${chem_arg} \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()} \\
        --nosecondary
    """
}

// Cell Ranger vdj (BCR)
process CELLRANGER_VDJ {
    tag "${sample_id}"
    label 'process_high'
    publishDir "${params.outdir}/cellranger_vdj/${sample_id}", mode: 'copy'

    container params.cellranger_container

    input:
    tuple val(sample_id), path(fastq_dir)
    path  vdj_reference

    output:
    tuple val(sample_id), path("${sample_id}/outs/filtered_contig_annotations.csv"),
        emit: contigs
    tuple val(sample_id), path("${sample_id}/outs/clonotypes.csv"),
        emit: clonotypes
    tuple val(sample_id), path("${sample_id}/outs/web_summary.html"),
        emit: web_summary
    tuple val(sample_id), path("${sample_id}/outs/metrics_summary.csv"),
        emit: metrics

    script:
    """
    cellranger vdj \\
        --id=${sample_id} \\
        --reference=${vdj_reference} \\
        --fastqs=${fastq_dir} \\
        --sample=${sample_id}_BCR \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()}
    """
}

// Helper: Build sample channels
def build_gex_channel() {
    // Get all samples WITHOUT 'BCR' in the name
    return Channel
        .fromPath("${params.fastq_dir}/*_S*_L001_R1_001.fastq.gz")
        .map { f ->
            // Extract sample name: Control1_S1_L001_R1_001.fastq.gz → Control1
            def matcher = (f.name =~ /^(.+?)_S\d+_L\d+_R\d+_\d+\.fastq\.gz$/)
            if (matcher) {
                def sample_id = matcher[0][1]
                if (!sample_id.contains('BCR')) {
                    return tuple(sample_id, f.parent)
                }
            }
            return null
        }
        .filter { it != null }
        .unique { it[0] }  // unique by sample_id
}

def build_bcr_channel() {
    // Get all samples WITH 'BCR' in the name
    return Channel
        .fromPath("${params.fastq_dir}/*_BCR_S*_L001_R1_001.fastq.gz")
        .map { f ->
            // Extract sample name: Control1_BCR_S9_L001_R1_001.fastq.gz → Control1
            def matcher = (f.name =~ /^(.+?)_BCR_S\d+_L\d+_R\d+_\d+\.fastq\.gz$/)
            if (matcher) {
                def sample_id = matcher[0][1]
                return tuple(sample_id, f.parent)
            }
            return null
        }
        .filter { it != null }
        .unique { it[0] }  // unique by sample_id
}

// Main Workflow
workflow {

    log.info """
    ╔══════════════════════════════════════════╗
    ║  Cell Ranger Dual Pipeline (GEX + BCR)   ║
    ╠══════════════════════════════════════════╣
    ║  fastq_dir    : ${params.fastq_dir}
    ║  transcriptome: ${params.transcriptome}
    ║  vdj_reference: ${params.vdj_reference}
    ║  chemistry    : ${params.chemistry}
    ║  outdir       : ${params.outdir}
    ╚══════════════════════════════════════════╝
    """.stripIndent()

    // Check inputs
    if (!params.fastq_dir) {
        error "Please provide --fastq_dir /path/to/fastqs"
    }
    if (!params.transcriptome) {
        error "Please provide --transcriptome /path/to/cellranger_ref"
    }
    if (!params.vdj_reference) {
        error "Please provide --vdj_reference /path/to/vdj_ref"
    }

    // Build channels
    gex_ch = build_gex_channel()
    bcr_ch = build_bcr_channel()

    // FastQC on all FASTQ files
    fastq_ch = Channel
        .fromPath("${params.fastq_dir}/*_R{1,2}_*.fastq.gz")
        .map { f ->
            def sample_base = f.name.replaceAll(/_S\d+.*/, '')
            tuple(sample_base, f)
        }
        .groupTuple()

    FASTQC(fastq_ch)
    MULTIQC(FASTQC.out.reports.collect())

    // Run pipelines
    transcriptome_ch = Channel.value(file(params.transcriptome, checkIfExists: true))
    vdj_ref_ch       = Channel.value(file(params.vdj_reference, checkIfExists: true))

    CELLRANGER_COUNT(gex_ch, transcriptome_ch)
    CELLRANGER_VDJ(bcr_ch, vdj_ref_ch)

    // Print completion
    CELLRANGER_COUNT.out.web_summary.view { sample_id, html ->
        "GEX ${sample_id} complete"
    }
    CELLRANGER_VDJ.out.web_summary.view { sample_id, html ->
        "BCR ${sample_id} complete"
    }
}

// Process resources
process {
    withLabel: 'process_low' {
        cpus   = 2
        memory = '8.GB'
        time   = '2.h'
    }
    withLabel: 'process_high' {
        cpus   = 32
        memory = '128.GB'
        time   = '48.h'
    }
    errorStrategy = { task.exitStatus in [130, 137, 140] ? 'retry' : 'finish' }
    maxRetries    = 2
}

// Execution profiles
profiles {
    
    slurm {
        process.executor = 'slurm'
        process.queue    = 'normal'
        singularity.enabled    = true
        singularity.autoMounts = true
        singularity.cacheDir   = params.container_cache
    }
    
    standard {
        process.executor = 'local'
        singularity.enabled    = true
        singularity.autoMounts = true
        singularity.cacheDir   = params.container_cache
    }
}