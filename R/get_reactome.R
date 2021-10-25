
library(AnnotationDbi)
library(org.Hs.eg.db)

# Download Reactome database
# NCBI Entrez returns 1:1 mapping to gene symbols
rea <- readLines('https://reactome.org/download/current/NCBI2Reactome_All_Levels.txt')
rea_l <- lapply(rea, function(i) strsplit(i, "\t")[[1]])
rea_df <- do.call(rbind, rea_l)
rea_df <- rea_df[rea_df[,6] == "Homo sapiens",]

symbols <- mapIds(org.Hs.eg.db, rea_df[,1], 'SYMBOL', 'ENTREZID')
rea_df <- cbind(symbols, rea_df)
rea_df <- rea_df[!is.na(rea_df[,1]),]

# genes to terms
genes_dict <- by(data = rea_df, INDICES = rea_df[,1], FUN = function(gene_df){
                 terms <- paste(gene_df[,3], collapse=',')
                 return(terms) })
genes_dict <- as.list(genes_dict)
genes2terms <- as.data.frame(do.call(rbind, genes_dict), stringsAsFactors=FALSE)
colnames(genes2terms) <- "Terms"

# terms to genes
terms_dict <- by(data = rea_df, INDICES = rea_df[,3], FUN = function(term_df){
                 genes <- paste(term_df[,1], collapse=',')
                 return(genes) })
terms_dict <- as.list(terms_dict)
terms2genes <- as.data.frame(do.call(rbind, terms_dict), stringsAsFactors=FALSE)
colnames(terms2genes) <- "Genes"

# terms descriptions
path2sum <- readLines('https://reactome.org/download/current/pathway2summation.txt')
path2sum_l <- lapply(path2sum, function(i) strsplit(i, "\t")[[1]])
path2sum_df <- do.call(rbind, path2sum_l)
k <- !duplicated(path2sum_df[,1])
path2sum_df <- path2sum_df[k,]
rownames(path2sum_df) <- path2sum_df[,1]
path2sum_df <- path2sum_df[,-1]
colnames(path2sum_df) <- c("Name", "Description")

k <- rownames(path2sum_df) %in% rownames(terms2genes)
path2sum_df <- path2sum_df[k,]

#########################################################################
# Write output
#########################################################################

dir_out <- "data/human/Reactome/"
dir.create(dir_out, showWarnings=FALSE, recursive=TRUE)

outfn <- paste0(dir_out, "term_description.txt")
write.table(path2sum_df, outfn, row.names=TRUE, col.names=NA, sep="\t", quote=TRUE)

outfn <- paste0(dir_out, "genes2terms.txt")
write.table(genes2terms, outfn, row.names=TRUE, col.names=NA, sep="\t", quote=FALSE)

outfn <- paste0(dir_out, "terms2genes.txt")
write.table(terms2genes, outfn, row.names=TRUE, col.names=NA, sep="\t", quote=FALSE)

out.l <- list("Source" = "Reactome",
              "Species" = "Human", 
              "genes2terms" = genes2terms, 
              "terms2genes" = terms2genes, 
              "term_description" = path2sum_df)

outfn <- paste0(dir_out,  "pathway.rds")
saveRDS(out.l, outfn)

# Save download date
date_df <- date()
write.table(date_df, paste0(dir_out, "download_date.txt"), 
	row.names=FALSE, col.names=FALSE, quote=TRUE)

