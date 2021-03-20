
library(KEGGREST)

#' parse genes from kegg list
parse_genes <- function(x){
    genes <- sapply(x, function(g) strsplit(g, "[,;]")[[1]][1])
    return(genes)
}

#' download all genes
#' org is organism, 'hsa' for human
all_genes <- function(org){
    require(KEGGREST)
    g <- keggList(org)
    g <- parse_genes(g)
    return(g)
}


get_pathway_genes <- function(org){
    require(KEGGREST)

    pathways <- keggList("pathway", org)

    pathway_list <- list()
    s <- seq(1, length(pathways), 10)
    s[length(s) + 1] <- length(pathways)+1

    for (i in 1:(length(s)-1)){
        Sys.sleep(0.1)
        query <- keggGet(names(pathways)[s[i]:(s[i+1] - 1)])
        for (q in query){
            store <- q[c("ENTRY", "NAME", "DESCRIPTION", "CLASS", "GENE")]
            if (length(store$GENE) <= 0) next
            store$GENE <- parse_genes(store$GENE[seq(2, length(store$GENE), 2)])
            pathway_list[[store$ENTRY]] <- store
        }
    }
    return(pathway_list)
}

org <- "hsa"

hsa_genes <- all_genes("hsa")
hsa_genes <- unique(hsa_genes)
pathway_list <- get_pathway_genes("hsa")

path_genes <- c()
for (p in pathway_list){
    path_genes <- c(path_genes, p$GENE)
}
path_genes <- unique(path_genes)

# genes to pathway
genes2pathway.l <- list()
for (p in pathway_list){
    pathway <- p$ENTRY
    for (gene in p$GENE){
        genes2pathway.l[[gene]] = c(genes2pathway.l[[gene]],pathway)
    }
}

genes2pathway.df <- data.frame("Gene" = names(genes2pathway.l), 
                               "Pathways" = sapply(genes2pathway.l, function(p){
                                                   paste(p, collapse=",") })
                               )

# pathway to genes
pathway2genes.df <- data.frame("Pathway" = names(pathway_list), 
                               "Genes" = sapply(pathway_list, function(p){
                                                paste(p$GENE, collapse=",") })
                               )

# pathway info
pathways.df <- data.frame("Pathway" = names(pathway_list), 
                          "Name" = sapply(pathway_list, function(p) p$NAME), 
                          "Description" = sapply(pathway_list, function(p) paste0(p$DESCRIPTION, collapse=';')))



# write output
dir_out <- "data/ref/KEGG/"
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

out.l <- list("pathway.list" = pathway_list, 
              "genes2pathway" = genes2pathway.df, 
              "pathway2genes" = pathway2genes.df, 
              "pathways" = pathways.df)
outfn <- paste0(dir_out,  "KEGG.rds")
saveRDS(out.l, outfn)

outfn <- paste0(dir_out,  "genes2terms.txt")
write.table(genes2pathway.df, outfn, row.names=FALSE, col.names=TRUE, 
            quote=FALSE, sep="\t")

outfn <- paste0(dir_out,  "terms2genes.txt")
write.table(pathway2genes.df, outfn, row.names=FALSE, col.names=TRUE, 
            quote=FALSE, sep="\t")

outfn <- paste0(dir_out,  "term_description.txt")
write.table(pathways.df, outfn, row.names=FALSE, col.names=TRUE, 
            quote=FALSE, sep="\t")

# outfn <- paste0(dir_out,  "genes.human.txt")
# writeLines(hsa_genes, outfn)

date_df <- date()
write.table(date_df, paste0(dir_out, "download_date.txt"), 
	row.names=FALSE, col.names=FALSE, quote=TRUE)


