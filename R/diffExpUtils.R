'%notin%' <- function(x,y)!('%in%'(x,y))
cbbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


write_DESeq2_results <- function(df, results.dir, prefix){
    df<- as.data.frame(df)
    df <- df[order(df$padj),]
    df$gene_name <- gene_annotations[rownames(df),]$gene_name

    df.sig <- subset(df, padj<0.05)
    df.sig.up <- subset(df.sig, log2FoldChange>0)
    df.sig.down <- subset(df.sig, log2FoldChange<0)
    write.table(df, file = file.path(results.dir,
                                     paste(prefix, 'tsv', sep='.')), sep = '\t')

    write.table(df.sig, file = file.path(results.dir,
                                         paste(prefix, 'sig', 'tsv', sep='.')), sep = '\t')
    write.table(df.sig.up,  file = file.path(results.dir,
                                             paste(prefix, 'sig', 'up', 'tsv', sep='.')), sep = '\t')
    write.table(df.sig.down,  file = file.path(results.dir,
                                               paste(prefix, 'sig', 'down', 'tsv', sep='.')), sep = '\t')
    return (df.sig)
}

plotHeatMap <- function(rlogdist, filename=NULL){
    sampleDists <- dist(t(assay(rlogdist)))
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- colnames(rlogdist)#paste(rlogdist$condition, colnames(rlogdist), sep="-")
    colnames(sampleDistMatrix) <- NULL
    colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
    pheatmap(sampleDistMatrix,
             clustering_distance_rows=sampleDists,
             cellwidth=10,
             cellheight=10,
             clustering_distance_cols=sampleDists,
             col=colors,)
    if (!is.null(filename)) {
        pheatmap(sampleDistMatrix,
                 cellwidth=20,
                 cellheight=20,
                 clustering_distance_rows=sampleDists,
                 clustering_distance_cols=sampleDists,
                 col=colors, filename=filename)
    }
}


readFcountsMatrix <- function(filepath, remove_txid_version=TRUE){
    counts.data <- read.table(filepath,
                              row.names=1,
                              stringsAsFactors = F,
                              header = T)

    if (remove_txid_version){
        rownames(counts.data) <- gsub('\\.[0-9]+', '', rownames(counts.data))
        }
    # Remove .bam or .sam from filenames
    colnames(counts.data) <- gsub("\\.[sb]am$", "", colnames(counts.data))
    colnames(counts.data) <- gsub("bams_unique\\.", "", colnames(counts.data))
    return (counts.data)

}

plotDESeq2PCA <- function(rld, intgroup=c("condition", "assay"), colorBy, shapeBy, colors){
    data <- plotPCA(rld, intgroup = intgroup, returnData=TRUE)
    percentVar <- round(100 * attr(data, "percentVar"))

    p<- ggplot(data, aes_string("PC1", "PC2", color=colorBy, shape=shapeBy, label = rownames(data))) +
        scale_color_manual(intgroup[1], values=colors) +
        geom_text_repel() +
        geom_point(size=3, alpha=0.3) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) +
        coord_fixed() + #ggtitle('PCA -- all samples') +
        theme(text = element_text(size=12))
    return (p)
}

doPvalueAdjustment <- function(results){
    hist(results$pvalue,  main = 'DESeq2 unadjusted p-values',
         xlab = 'Unadjusted p-values')
    results <- results[ !is.na(results$padj), ]
    results <- results[ !is.na(results$pvalue), ]
    results <- results[, -which(names(results) == 'padj')]
    resultsFDR <- fdrtool(results$stat,
                          statistic= 'normal',
                          plot = T)
    results[,'padj']  <- p.adjust(resultsFDR$pval,
                                  method = 'BH')
    results[,'padj.renull'] <-   resultsFDR$pval
    hist(resultsFDR$pval,
         main = 'DESeq2 corrected p-values | Empirical null',
         xlab = 'Corrected p-values')
    return (results)
}
