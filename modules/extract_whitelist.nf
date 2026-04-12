process EXTRACT_WHITELIST {

    tag "$meta.id"
    label 'process_low'

    publishDir "${params.outdir}/${params.genome}/cellranger/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(matrix_dir)

    output:
    tuple val(meta), path("${meta.id}_whitelist.txt"), emit: whitelist

    script:
    """
    zcat ${matrix_dir}/barcodes.tsv.gz | sed 's/-1\$//' > ${meta.id}_whitelist.txt

    NBARCODES=\$(wc -l < ${meta.id}_whitelist.txt)
    if [ "\$NBARCODES" -eq 0 ]; then
        echo "ERROR: Whitelist for ${meta.id} is empty." >&2
        exit 1
    fi
    echo "[EXTRACT_WHITELIST] ${meta.id}: \$NBARCODES barcodes written."
    """
}