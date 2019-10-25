#' Convert Entrez ID to Gene symbol (Human)
#'
#' @param entrezID Entrez ID (Human),
#' @return Gene Symbol
#'
#'
EntrezToSymbolHuman <- function(entrezID){
    bitr(entrezID, fromType = 'ENTREZID', toType = 'SYMBOL', OrgDb = org.Hs.eg.db)$SYMBOL
}

#' Write clusterProfile KEGG data to csv
#'
#' @param results ekg clusterProfiler enrichKEGG object.
#' @param csvout Path to output csv,
#' @param converter Function to map ENtrez IDs to Gene Symbol (Default=EntrezToSymbolHuman)
#'
#'
clusterProfilerWriteKEGG <- function(ekg, csvout, converter = EntrezToSymbolHuman){
    ekg <- as.data.frame(ekg)
    ekg.geneID <- strsplit(as.character(ekg$geneID), '/', fixed=TRUE)
    ekg.genes <- lapply(ekg.geneID, converter)
    ekg$geneID <- ekg.genes
    ekg$geneID <- vapply(ekg$geneID , paste, collapse = ", ", character(1L))
    write.csv(ekg, csvout)
}

#' Write clusterProfile GO data to csv
#'
#' @param results ekg clusterProfiler enrichGO object.
#' @param csvout Path to output csv,
#' @param converter Function to map ENtrez IDs to Gene Symbol (Default=EntrezToSymbolHuman)
#'
#'
clusterProfilerWriteGO <- function(ego, csvout, converter = EntrezToSymbolHuman){
    ego <- as.data.frame(ego)
    ego.geneID <- strsplit(as.character(ego$geneID), '/', fixed=TRUE)
    ego.genes <- lapply(ego.geneID, converter)
    ego$geneID <- ego.genes
    ego$geneID <- vapply(ego$geneID , paste, collapse = ", ", character(1L))
    write.csv(ego, csvout)
}



doTopGoAnalysis <- function(geneID2GO, interestingGenes, ontology="BP", universalGeneList="allInGO", topNodes=15){
    if (universalGeneList == "allInGO"){
        universalGeneList <- names(geneID2GO)
    }

    geneList <- factor(as.integer(universalGeneList %in% interestingGenes))
    names(geneList) <- universalGeneList
    GOdata <- new("topGOdata",
                  ontology = ontology,
                  allGenes = geneList,
                  annot = annFUN.gene2GO,
                  gene2GO = geneID2GO)
    resultFisher <- runTest(GOdata,
                            algorithm = "classic",
                            statistic = "fisher")
    resultKS <- runTest(GOdata,
                        algorithm = "classic",
                        statistic = "ks")
    resultKS.elim <- runTest(GOdata,
                             algorithm = "elim",
                             statistic = "ks")
    allRes <- GenTable(GOdata,
                       classicFisher = resultFisher,
                       classicKS = resultKS,
                       elimKS = resultKS.elim,
                       orderBy = "elimKS",
                       ranksOf = "classicFisher",
                       topNodes = topNodes)
    return (as.data.frame(allRes))


}


barplotGOSeq <- function(df, showCategory=15){
    df <- df[with(df, order(ratio, padj, decreasing = c(TRUE, FALSE))),]
    df <- head(df, n=showCategory)
    breaks <- round( c(0, 1/4, 2/4, 3/4, 1) * max(df[['ratio']]) , 2)
    p_plot <- ggplot(df, aes_string(x="term", y="ratio", fill="padj")) +
        geom_col() +
        scale_y_continuous(expand=c(0, 0), breaks=breaks, limits=c(0, max(df[["ratio"]]+0.05))) +
        scale_x_discrete(name='GO term') +
        scale_fill_continuous(low="#00dbde", high="#FFF94C") +
        theme(text=ggplot2::element_text(size=9)) +
        coord_flip() +
        theme_bw(base_size=9)
    return(p_plot)
}

dotplotGOSeq <- function(df, showCategory=15){
    df <- df[with(df, order(ratio, padj, decreasing = c(TRUE, FALSE))),]
    df <- head(df, n=showCategory)
    d_plot <- ggplot(df, aes_string(x="term",
                                    y="ratio",
                                    colour="padj",
                                    size="numDEInCat")) +
        geom_point() +
        scale_color_gradient(low="#00dbde",
                             high="#FFF94C") +
        coord_flip() +
        theme_bw(base_size=9)
    return(d_plot)
}
