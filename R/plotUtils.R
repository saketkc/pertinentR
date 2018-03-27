#' Draw volcano plot using DESeq2 results object.
#'
#' @param results DESeq2 results obj.
#' @param numerator Name of numerator  of contrast
#' @param denominator Name of denominator  of contrast
#' @param pValCutoff Cutoff to draw the lines
#' @return plotObject
#'
#'
#'
library()
cbbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
plotVolcanoDESeq <- function(results, numerator, denominator, pValCutoff=1){
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
