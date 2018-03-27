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
