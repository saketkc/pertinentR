cbbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#' Draw volcano plot using DESeq2 results object.
#'
#' @param results DESeq2 results obj.
#' @param numerator Name of numerator  of contrast
#' @param denominator Name of denominator  of contrast
#' @param pValCutoff Cutoff to draw the lines
#' @return plotObject
#'
#'
DESeqPlotVolcano <- function(results, numerator, denominator, pValCutoff=1){
    data <- data.frame(gene = row.names(results),
                       pvalue = -log10(results$padj),
                       lfc = results$log2FoldChange)

    # Remove rows that have NA values
    data <- na.omit(data)
    data <- data %>%  mutate(color = ifelse(data$lfc > 0 & data$pvalue > pValCutoff,
                                            yes = numerator,
                                            no = ifelse(data$lfc < 0 & data$pvalue > pValCutoff,
                                                        yes = denominator,
                                                        no = 'none')))
    colors <- c()
    colors[numerator] <- cbbPalette[6]
    colors[denominator] <- cbbPalette[7]
    colors[['none']] <- cbbPalette[1]
    labels <- c()
    labels[numerator] <- numerator
    labels[denominator] <- denominator
    labels['none'] <- 'n.s.'
    xlabel <- bquote(log[2](.(numerator)/.(denominator)))

    colored <- ggplot(data, aes(x = lfc, y = pvalue)) +
        geom_point(aes(color = factor(color)), size = 0.5, alpha = 0.5, na.rm = T) + # add gene points
        theme_bw(base_size = 16) + # clean up theme
        theme(legend.position = 'none', plot.title = element_text(hjust = 0.5)) + # remove legend
        ylab(expression(-log[10]('adjusted p-value'))) +
        geom_vline(xintercept = 0, colour = 'black', linetype='dashed') +
        geom_hline(yintercept = pValCutoff, colour = 'black', linetype='dashed')  +
        xlab(xlabel) +
        theme(aspect.ratio=1) +
        coord_fixed() +
        scale_color_manual(values = colors) +
        labs(color = "")

    return (colored)
}

#' Draw sample to sample heatmap using DESeq2 object.
#'
#' @param dds.vst DESeq2 post VST obj .
#' @param filename Filename to save to
#'
#'
DESeqPlotSampleDistmap <- function (dds.vst, filename=NULL){
    sampleDists <- dist(t(assay(dds.vst)))
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- colnames(dds.vst)
    colnames(sampleDistMatrix) <- NULL
    colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
    if (is.null(filename)){
        pheatmap(sampleDistMatrix,
                 clustering_distance_rows = sampleDists,
                 clustering_distance_cols = sampleDists,
                 col = colors,
                 title = title,
                 fontsize = 2)
    }
    else {
        pheatmap(sampleDistMatrix,
                 clustering_distance_rows = sampleDists,
                 clustering_distance_cols = sampleDists,
                 col=colors,
                 title=title,
                 fontsize=2,
                 filename=filename)
    }
}

#' Draw PCA using DESeq2 object.
#'
#' @param dds.vst DESeq2 post VST obj .
#' @param intrgroup Features to plot (intgroup[1] - color, intgroup[2]-shape)
#'
#'

DESeqPlotPCA <- function(dds.vst, intgroup, title){
    data <- plotPCA(dds.vst, intgroup = intgroup, returnData=TRUE)
    percentVar <- round(100 * attr(data, "percentVar"))
    p <- ggplot(data, aes(PC1, PC2, color=bquote(intgroup[1]), shape = intgroup[2])) +
        geom_point(size=5, alpha=0.7) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) +
        coord_fixed() + ggtitle(title) +
        theme(text = element_text(size=12))
    return (p)
}

#' Draw sample to sample heatmap using DESeq2 object.
#'
#' @param dds.vst DESeq2 post VST obj .
#' @param col_df Design dataframe
#' @param col_labels Label names for columns
#' @param title Plot title
#' @param filename Path to file
#'
DESeqPlotHeatmap <- function(dds, col_df, col_labels, title, filename){
    #df <- as.data.frame(colData(dds)[, c("condition","age", "sex", "estimated_batch")])
    #col_labels <- paste(design.info$condition, design.info$sampleID)
    pheatmap(assay(dds),
             cluster_rows = FALSE,
             show_rownames = FALSE,
             cluster_cols = TRUE,
             annotation_col = col_df,
             fontsize = 4,
             title = title,
             labels_col = col_labels,
             filename = filename)
}
