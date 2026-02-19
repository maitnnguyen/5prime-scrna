// ============================================================
//  modules/multiqc.nf
//  Aggregate all QC reports into a single HTML report
// ============================================================

process MULTIQC {

    label 'process_low'
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path logs   // All log/QC files collected from every module

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data/",       emit: data

    script:
    def config = params.multiqc_config ? "--config ${params.multiqc_config}" : ""
    def title  = params.multiqc_title  ? "--title \"${params.multiqc_title}\"" : ""
    """
    multiqc \\
        ${config} \\
        ${title} \\
        --force \\
        --outdir . \\
        .
    """
}
