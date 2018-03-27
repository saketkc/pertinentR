
clusterProfilergetGeneName <- function(x){
    bitr(x, fromType = 'ENTREZID', toType = 'SYMBOL', OrgDb = org.Hs.eg.db)$SYMBOL
}

clusterProfilerWriteKEGG <- function(ekg, csvout){
    ekg <- as.data.frame(ekg)
    ekg.geneID <- strsplit(as.character(ekg$geneID), '/', fixed=TRUE)
    ekg.genes <- lapply(ekg.geneID, get_name)
    ekg$geneID <- ekg.genes
    ekg$geneID <- vapply(ekg$geneID , paste, collapse = ", ", character(1L))
    write.csv(ekg, csvout)
}

clusterProfilerWriteGO <- function(ego, csvout){
    ego <- as.data.frame(ego)
    ego.geneID <- strsplit(as.character(ego$geneID), '/', fixed=TRUE)
    ego.genes <- lapply(ego.geneID, get_name)
    ego$geneID <- ego.genes
    ego$geneID <- vapply(ego$geneID , paste, collapse = ", ", character(1L))
    write.csv(ego, csvout)
}
