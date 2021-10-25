
library(KEGGREST)

#########################################################################
# Functions
#########################################################################

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

#########################################################################
# Download data
#########################################################################

org <- "hsa"

hsa_genes <- all_genes("hsa")
hsa_genes <- unique(hsa_genes)
pathway_list <- get_pathway_genes("hsa")

#########################################################################
# Write output
#########################################################################

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

genes2pathway.df <- data.frame("Terms" = sapply(genes2pathway.l, function(p){
                                                   paste(p, collapse=",") }), 
                               stringsAsFactors=FALSE)
rownames(genes2pathway.df) <- names(genes2pathway.l)

# pathway to genes
pathway2genes.df <- data.frame("Genes" = sapply(pathway_list, function(p){
                                                paste(p$GENE, collapse=",") }),
                               stringsAsFactors=FALSE)
rownames(pathway2genes.df) <- names(pathway_list)

# pathway info
pathways.df <- data.frame("Name" = sapply(pathway_list, function(p) p$NAME), 
                          "Description" = sapply(pathway_list, function(p) paste0(p$DESCRIPTION, collapse=';')), 
                          stringsAsFactors=FALSE)
rownames(pathways.df) <- names(pathway_list)

# get category
pathways.df[,"Class"] <- sapply(pathway_list[rownames(pathways.df)], 
                                   function(x){ x[["CLASS"]] })

# write output
dir_out <- "data/human/KEGG/"
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

out.l <- list("Source" = "KEGG", 
              "Species" = "Human", 
              "genes2terms" = genes2pathway.df, 
              "terms2genes" = pathway2genes.df, 
              "term_description" = pathways.df)
outfn <- paste0(dir_out,  "pathway.rds")
saveRDS(out.l, outfn)

outfn <- paste0(dir_out,  "genes2terms.txt")
write.table(genes2pathway.df, outfn, row.names=TRUE, col.names=NA, 
            quote=FALSE, sep="\t")

outfn <- paste0(dir_out,  "terms2genes.txt")
write.table(pathway2genes.df, outfn, row.names=TRUE, col.names=NA, 
            quote=FALSE, sep="\t")

outfn <- paste0(dir_out,  "term_description.txt")
write.table(pathways.df, outfn, row.names=TRUE, col.names=NA, 
            quote=TRUE, sep="\t")

# outfn <- paste0(dir_out,  "genes.human.txt")
# writeLines(hsa_genes, outfn)

date_df <- date()
write.table(date_df, paste0(dir_out, "download_date.txt"), 
	row.names=FALSE, col.names=FALSE, quote=TRUE)


