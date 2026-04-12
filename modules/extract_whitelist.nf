process EXTRACT_WHITELIST {
    tag "${sample_id}"
    label 'process_low'

    input:
    tuple val(sample_id), path(matrix_dir) // from CELLRANGER_COUNT.out.filtered_matrix

    output:
    tuple val(sample_id), path("${sample_id}_whitelist.txt"), emit: whitelist

    script:
    """
    # Cell Ranger stores barcodes in filtered_feature_bc_matrix/barcodes.tsv.gz
    zcat ${matrix_dir}/barcodes.tsv.gz | sed 's/-1//' > ${sample_id}_whitelist.txt
    """
}