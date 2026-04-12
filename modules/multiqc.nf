// ============================================================
//  modules/multiqc.nf
//  Aggregate all QC reports into a single HTML report
// ============================================================

process MULTIQC {
    label 'process_low'
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    container params.multiqc_container

    input:
    path 'input_files??/*' // Using a glob pattern helps MultiQC find nested files

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data/",       emit: data

    script:
    """
    multiqc . \\
        --filename multiqc_report.html \\
        --force \\
        --interactive
    """
}