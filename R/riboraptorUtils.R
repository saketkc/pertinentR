#' Extend any vector to specified length filled with NAs
#'
#' @param vector Input vector.
#' @param newLength Final length of extended vecotr
#' @return Extended vector
#'
extendLength <- function(vector, newLength){
    origLength <- length(vector)
    vector <- rep(vector, length.out = newLength)
    vector[origLength:newLength] <- NA
    return(vector)
}

#' Parse count column of dataframe saved by riboraptor
#'
#' @param col Count column,
#' @param newLength Final length of extended vecotr
#' @return Extended vector post columns parsing
#'
riboSummaryParseCountColumn <- function(col, newLength){
    col <- gsub(']', ')', col)
    col <- gsub('\\[', 'c(', col)
    col <- eval(parse(text = col))
    return(extendLength(col, newLength))
}

#' Extend any vector to specified length filled with NAs
#'
#' @param df Data frame as read from riboSummary Counts file
#' @return Counts Matrix
#'
#'
riboSummaryCollapseDF <- function(df){
    df <- arrange(df,  -mean, -length)
    maxLength <- max(df$length)
    countsMatrix <- t(simplify2array((lapply(df$count, riboSummaryParseCountColumn, newLength = maxLength))))
    return(countsMatrix)
}
