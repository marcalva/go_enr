

#' Return keggGet output for pathway ID
get_kegg_path <- function(id){
    require(KEGGREST)
    g <- keggGet(id)
    return(g)
}

#' Get genes from a kegg database entry list
kegg_genes <- function(dbl){
    g <- dbl[['GENE']]
    gdf <- matrix(g, ncol = 2, byrow = TRUE)
    genes <- sapply(gdf[,2], function(x) strsplit(x, ';')[[1]][1])
    names(genes) <- NULL
    return(genes)
}

