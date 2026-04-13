// ============================================================
//  modules/ctss_counts_bigwig.nf
//  Count CTSS reads at FANTOM5 promoters/enhancers
//  Generate CPM-normalised BigWig for IGV / UCSC visualisation
//  Adapted from: STARsolo_Counts_CTSS_bed_CPM_bigwig_240119.sh
// ============================================================

process CTSS_COUNTS_BIGWIG {
    tag "$meta.id"
    label 'process_medium'
    
    // Organizes by genome (T2T/hg38)
    publishDir "${params.outdir}/${params.genome}/reaptec/counts_and_bw", mode: 'copy'

    input:
    tuple val(meta), path(ctss_bed)  // From REAPTEC_CTSS_BED
    path chrom_sizes                 // Sorted chrom.sizes file
    path promoter_ref                // FANTOM5 promoter bed
    path enhancer_ref                // FANTOM-NET enhancer bed

    output:
    tuple val(meta), path("${meta.id}.CTSS.CPM.fwd.bw"), emit: bw_fwd
    tuple val(meta), path("${meta.id}.CTSS.CPM.rev.bw"), emit: bw_rev
    path("${meta.id}.promoter.fwd.rev.txt"), optional: true, emit: prom_counts
    path("${meta.id}.enhancer.*.txt"),       optional: true, emit: enh_counts

    script:
    """
    # load bedtools module
    module load BEDTools/2.31.1-GCC-14.3.0
    module load OpenSSL/3.6
    
    # 1. Calculate Total CTSS sum for CPM (Column 5 in your CTSS.bed)
    SUM_TOTAL=\$(awk 'BEGIN{sum=0}{sum=sum+\$5}END{print sum}' ${ctss_bed})

    # 2. Collapse CTSS to avoid overlapping regions (Murakawa Script Step)
    awk 'BEGIN {OFS="\\t"} {print \$1, \$2, \$3, ".", \$5, \$6}' ${ctss_bed} | \\
        bedtools groupby -g 1,2,3,4,6 -c 5 -o sum | \\
        awk 'BEGIN{OFS="\\t"}{print \$1, \$2, \$3, \$4, \$6, \$5}' > collapsed.bed

    # 3. Create CPM BigWigs
    # Forward (+)
    awk -v sum=\$SUM_TOTAL 'BEGIN{OFS="\\t"} \$5 > 0 && \$6=="+" {printf("%s\\t%i\\t%i\\t%1.2f\\n", \$1,\$2,\$3, 1e6 * \$5 / sum)}' collapsed.bed | \\
        sort -k1,1 -k2,2n > fwd.bg
    bedGraphToBigWig fwd.bg ${chrom_sizes} ${meta.id}.CTSS.CPM.fwd.bw

    # Reverse (-)
    awk -v sum=\$SUM_TOTAL 'BEGIN{OFS="\\t"} \$5 > 0 && \$6=="-" {printf("%s\\t%i\\t%i\\t%1.2f\\n", \$1,\$2,\$3, 1e6 * \$5 / sum)}' collapsed.bed | \\
        sort -k1,1 -k2,2n > rev.bg
    bedGraphToBigWig rev.bg ${chrom_sizes} ${meta.id}.CTSS.CPM.rev.bw

    # 4. FANTOM Quantification (Only if files are provided and non-empty)
    # This logic matches the bash script exactly
    if [[ -f "${promoter_ref}" && -s "${promoter_ref}" ]]; then
        # Prepare refs with 0 in score col
        awk 'BEGIN{OFS="\\t"}{print \$1, \$2, \$3, \$4, "0", \$6}' ${promoter_ref} > Promoter.tmp.bed
        awk 'BEGIN{OFS="\\t"}{print \$1, \$2, \$3, \$4, "0", \$6}' ${enhancer_ref} > Enhancer.tmp.bed

        # Promoter Counts
        cat ${ctss_bed} Promoter.tmp.bed | sort -k 1,1 -k 2,2n | \\
            bedtools merge -d -1 -s -c 5,6 -o sum,distinct -i stdin | \\
            awk 'BEGIN{OFS="\\t"}{print \$1, \$2, \$3, ".", \$4, \$5}' | \\
            bedtools intersect -wa -wb -s -a Promoter.tmp.bed -b stdin | \\
            cut -f 4,11 | sort > ${meta.id}.promoter.fwd.rev.txt

        # Enhancer Counts (Forward)
        awk '\$6=="+"' ${ctss_bed} | cat - Enhancer.tmp.bed | sort -k 1,1 -k 2,2n | \\
            bedtools merge -d -1 -c 5 -o sum -i stdin | \\
            awk 'BEGIN{OFS="\\t"}{print \$1, \$2, \$3, ".", \$4, "+"}' | \\
            bedtools intersect -wa -wb -a Enhancer.tmp.bed -b stdin | \\
            cut -f 4,11 | sort > ${meta.id}.enhancer.fwd.txt

        # Enhancer Counts (Reverse)
        awk '\$6=="-"' ${ctss_bed} | cat - Enhancer.tmp.bed | sort -k 1,1 -k 2,2n | \\
            bedtools merge -d -1 -c 5 -o sum -i stdin | \\
            awk 'BEGIN{OFS="\\t"}{print \$1, \$2, \$3, ".", \$4, "-"}' | \\
            bedtools intersect -wa -wb -a Enhancer.tmp.bed -b stdin | \\
            cut -f 4,11 | sort > ${meta.id}.enhancer.rev.txt
    fi
    """
}