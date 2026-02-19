// ============================================================
//  modules/seurat_clustering.nf
//  B Cell subtype clustering + enhancer matrix creation
// ============================================================

process SEURAT_CLUSTERING {

    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}/seurat/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(gene_matrix)       // STARsolo filtered/ dir
    tuple val(meta), path(barcodes)           // Filtered barcodes
    tuple val(meta), path(enhancers)          // btcEnh BED
    tuple val(meta), path(ctss_bed)           // CTSS BED for enhancer matrix

    output:
    tuple val(meta), path("${meta.id}.seurat.rds"),             emit: seurat_rds
    tuple val(meta), path("${meta.id}.cluster_barcodes/"),      emit: cluster_assignments
    tuple val(meta), path("${meta.id}.umap.pdf"),               emit: umap
    tuple val(meta), path("${meta.id}.seurat.log"),             emit: log

    script:
    """
    #!/usr/bin/env Rscript
    # ── Seurat B cell clustering pipeline ─────────────────────────

    suppressPackageStartupMessages({
        library(Seurat)
        library(dplyr)
        library(ggplot2)
    })

    log_file <- "${meta.id}.seurat.log"
    cat("[Seurat] Starting B cell clustering for: ${meta.id}\\n", file=log_file)

    # ── 1. Load STARsolo gene expression matrix ────────────────────
    cat("[Seurat] Loading gene expression matrix...\\n", file=log_file, append=TRUE)
    counts <- Read10X("${gene_matrix}")

    obj <- CreateSeuratObject(
        counts   = counts,
        min.cells     = ${params.min_cells},
        min.features  = ${params.min_features},
        project  = "${meta.id}"
    )

    cat(sprintf("[Seurat] Cells after filtering: %d\\n", ncol(obj)), file=log_file, append=TRUE)

    # ── 2. QC: mitochondrial gene filtering ────────────────────────
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern="^MT-")
    obj <- subset(obj, subset = percent.mt < ${params.max_mt_pct})

    cat(sprintf("[Seurat] Cells after MT filter (<${params.max_mt_pct}%%): %d\\n", ncol(obj)),
        file=log_file, append=TRUE)

    # ── 3. Normalization + variable features ───────────────────────
    obj <- NormalizeData(obj)
    obj <- FindVariableFeatures(obj, nfeatures=3000)
    obj <- ScaleData(obj)
    obj <- RunPCA(obj, npcs=${params.n_pcs})

    # ── 4. Clustering ──────────────────────────────────────────────
    obj <- FindNeighbors(obj, dims=1:${params.n_pcs})
    obj <- FindClusters(obj, resolution=${params.resolution})
    obj <- RunUMAP(obj, dims=1:${params.n_pcs})

    # ── 5. B cell subtype annotation ──────────────────────────────
    # Key B cell markers:
    #   Transitional: CD24hi, CD38hi, SELL
    #   Naïve:        CD27-, IGHD+, FCER2 (CD23)
    #   GC B cells:   BCL6+, AICDA+, CD38hi, CXCR4+
    #   Memory:       CD27+, IGHG/IGHA (switched) or IGHM (unswitched)
    #   Plasmablast:  PRDM1+, XBP1+, SDC1 (CD138+)
    bcell_markers <- list(
        "Transitional_B"  = c("CD24", "CD38", "SELL"),
        "Naive_B"         = c("FCER2", "IGHD", "TCL1A"),
        "GC_B"            = c("BCL6", "AICDA", "CXCR4", "MEF2B"),
        "Switched_Memory" = c("CD27", "IGHG1", "IGHG2", "IGHA1"),
        "Unswitched_Mem"  = c("CD27", "IGHM", "IGHD"),
        "Plasmablast"     = c("PRDM1", "XBP1", "IRF4", "SDC1")
    )

    # Score each cluster for each B cell subtype
    obj <- AddModuleScore(obj, features=bcell_markers, name="Bcell_score")

    cat("[Seurat] B cell subtype scoring complete.\\n", file=log_file, append=TRUE)

    # ── 6. UMAP plots ──────────────────────────────────────────────
    pdf("${meta.id}.umap.pdf", width=14, height=10)

    p1 <- DimPlot(obj, reduction="umap", group.by="seurat_clusters",
                  label=TRUE, pt.size=0.3) + ggtitle("${meta.id} - Seurat Clusters")
    print(p1)

    # B cell marker feature plots
    key_genes <- c("MS4A1","CD27","IGHD","BCL6","AICDA","PRDM1","XBP1","CD38")
    key_genes <- key_genes[key_genes %in% rownames(obj)]
    p2 <- FeaturePlot(obj, features=key_genes, ncol=4, pt.size=0.1)
    print(p2)

    dev.off()

    # ── 7. Export per-cluster barcode lists ───────────────────────
    dir.create("${meta.id}.cluster_barcodes", showWarnings=FALSE)

    clusters <- unique(obj\$seurat_clusters)
    for (cl in clusters) {
        cl_cells <- WhichCells(obj, idents=cl)
        writeLines(cl_cells,
            paste0("${meta.id}.cluster_barcodes/cluster_", cl, "_barcodes.txt"))
    }

    cat(sprintf("[Seurat] Exported barcodes for %d clusters.\\n", length(clusters)),
        file=log_file, append=TRUE)

    # ── 8. Save Seurat object ──────────────────────────────────────
    saveRDS(obj, "${meta.id}.seurat.rds")
    cat("[Seurat] Seurat object saved.\\n", file=log_file, append=TRUE)
    """
}


// ============================================================
//  modules/pseudobulk_enhancer.nf
//  Pseudobulk enhancer differential analysis (CVID vs Control)
// ============================================================

process PSEUDOBULK_ENHANCER {

    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}/pseudobulk/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(enhancers)
    tuple val(meta), path(cluster_dir)
    tuple val(meta), path(ctss_bed)
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.enhancer_counts.tsv.gz"),  emit: count_matrix
    tuple val(meta), path("${meta.id}.diff_enhancers.tsv.gz"),   emit: diff_results
    tuple val(meta), path("${meta.id}.pseudobulk.log"),           emit: log

    script:
    """
    #!/bin/bash
    set -euo pipefail

    echo "[Pseudobulk] Building per-cluster enhancer count matrix: ${meta.id}" > ${meta.id}.pseudobulk.log

    # ── For each B cell cluster: count CTSS reads at enhancers ────
    header="enhancer_id\\tchr\\tstart\\tend"

    for BARCODE_FILE in ${cluster_dir}/*.txt; do
        CLUSTER_NAME=\$(basename \$BARCODE_FILE .txt)

        # Filter CTSS to this cluster's barcodes
        zcat ${ctss_bed} | \\
            awk -v barcodes="\$BARCODE_FILE" '
            BEGIN { while ((getline line < barcodes) > 0) bc[line]=1 }
            { if (bc[\$4]) print }
            ' > tmp_cluster_ctss.bed

        # Count CTSS reads overlapping each enhancer
        bedtools coverage \\
            -a ${enhancers} \\
            -b tmp_cluster_ctss.bed \\
            -counts \\
            -s | \\
            awk '{print \$4"\\t"\$7}' > tmp_cluster_counts.txt

        header="\${header}\\t\${CLUSTER_NAME}"
        echo "[Pseudobulk] Cluster \$CLUSTER_NAME: \$(wc -l < tmp_cluster_ctss.bed) CTSS reads" >> ${meta.id}.pseudobulk.log
    done

    # ── Build count matrix ─────────────────────────────────────────
    # Add enhancer coordinates
    awk 'BEGIN{OFS="\\t"}{print \$4, \$1, \$2, \$3}' ${enhancers} > enh_coords.tsv

    # Merge all cluster counts
    paste enh_coords.tsv tmp_*_counts.txt | \\
        gzip > ${meta.id}.enhancer_counts.tsv.gz

    echo "[Pseudobulk] Enhancer count matrix built." >> ${meta.id}.pseudobulk.log

    # ── DESeq2 differential analysis (R) ──────────────────────────
    Rscript - <<'REOF'
    suppressPackageStartupMessages({
        library(DESeq2)
        library(dplyr)
        library(readr)
    })

    cat("[DESeq2] Running differential enhancer analysis...\\n")

    counts_mat <- read_tsv("${meta.id}.enhancer_counts.tsv.gz",
                           col_names=c("enhancer_id","chr","start","end"))

    # Extract numeric count columns
    count_data <- counts_mat[, -(1:4)]
    rownames(count_data) <- counts_mat\$enhancer_id

    # Build sample metadata from column names
    # Condition encoded in sample name: e.g. CVID_D1_cluster_0, HC_D2_cluster_2
    col_meta <- data.frame(
        sample    = colnames(count_data),
        condition = ifelse(grepl("CVID|cvid|patient", colnames(count_data)), "CVID", "Control"),
        row.names = colnames(count_data)
    )

    if (length(unique(col_meta\$condition)) < 2) {
        cat("[DESeq2] Only one condition found - skipping differential analysis.\\n")
        write_tsv(data.frame(note="Single condition - DESeq2 skipped"),
                  "${meta.id}.diff_enhancers.tsv.gz")
        quit(status=0)
    }

    # Run DESeq2
    dds <- DESeqDataSetFromMatrix(
        countData = round(as.matrix(count_data)),
        colData   = col_meta,
        design    = ~ condition
    )

    # Filter low-count enhancers
    dds <- dds[rowSums(counts(dds)) >= 10, ]

    dds <- DESeq(dds)
    res <- results(dds, contrast=c("condition","CVID","Control"),
                   alpha=0.05, lfcThreshold=0.5)

    # Add enhancer coordinates
    res_df <- as.data.frame(res)
    res_df\$enhancer_id <- rownames(res_df)
    res_df <- merge(res_df,
                    counts_mat[,c("enhancer_id","chr","start","end")],
                    by="enhancer_id")
    res_df <- res_df[order(res_df\$padj), ]

    # Save results
    write_tsv(res_df, "${meta.id}.diff_enhancers.tsv.gz")

    sig <- sum(res_df\$padj < 0.05 & abs(res_df\$log2FoldChange) > 0.5, na.rm=TRUE)
    cat(sprintf("[DESeq2] Significant differential enhancers (padj<0.05, |LFC|>0.5): %d\\n", sig))
    cat(sprintf("[DESeq2] CVID-gained enhancers: %d\\n",
        sum(res_df\$padj < 0.05 & res_df\$log2FoldChange > 0.5, na.rm=TRUE)))
    cat(sprintf("[DESeq2] CVID-lost enhancers: %d\\n",
        sum(res_df\$padj < 0.05 & res_df\$log2FoldChange < -0.5, na.rm=TRUE)))

REOF

    echo "[Pseudobulk] Analysis complete for: ${meta.id}" >> ${meta.id}.pseudobulk.log
    """
}
