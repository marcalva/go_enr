

#########################################################################
# Download and parse gene ontology terms, definitions, descriptions
#########################################################################

geo <- readLines("http://snapshot.geneontology.org/ontology/go-basic.obo")
# geo is a list, each element is a line of the file go-basic.obo

geo_dict <- list()
i <- 1
while (i <= length(geo)){
	if (geo[i] == "[Term]"){
		i <- i+1
		id <- sub("id: ", "", geo[i])
		i <- i+1
		name <- sub("name: ", "", geo[i])
		i <- i+1
		namespace <- sub("namespace: ", "", geo[i])
        i <- i+1
        def <- sub("def: ", "",  geo[i])
		geo_dict[[id]] <- c(name, namespace, def)
	}
	i <- i+1
}

geo_df <- as.data.frame(do.call(rbind, geo_dict))
colnames(geo_df) <- c("Name", "NameSpace", "Description")
desc_nq <- gsub("\"", "", geo_df[,"Description"])
geo_df[,"Description"] <- desc_nq
geo_df <- geo_df[,c("Name", "Description", "NameSpace")]

#########################################################################
# Genes to terms
#########################################################################

con <- gzcon(url("http://geneontology.org/gene-associations/goa_human.gaf.gz"))
txt <- readLines(con)
goa <- read.table(textConnection(txt), comment.char="!", header=FALSE, fill=TRUE, sep="\t",
				 stringsAsFactors=FALSE, quote="")

genes <- unique(goa[,3])
genes_dict <- by(data = goa, INDICES = factor(goa[,3]), FUN = function(gene_df){
                 uniprots = paste(unique(gene_df[,2]), collapse=",")
                 descr = paste(unique(gene_df[,10]), collapse=",")
                 go_terms = paste(unique(gene_df[,5]), collapse=",")
                 c(go_terms, descr, uniprots)})
genes_dict <- as.list(genes_dict)

genes_df <- as.data.frame(do.call(rbind, genes_dict), stringsAsFactors=FALSE)
colnames(genes_df) = c("Terms", "Description", "UniProtID")

genes2terms <- genes_df[,"Terms",drop=FALSE]
genes_descr <- genes_df[,c("Description", "UniProtID"),drop=FALSE]

#########################################################################
# Terms to genes
#########################################################################

terms <- unique(goa[,5])
terms_dict <- by(data = goa, INDICES = factor(goa[,5]), FUN = function(term_df){
                 genes <- paste(term_df[,3], collapse=",") })
terms_dict <- as.list(terms_dict)
                 
terms2genes <- as.data.frame(do.call(rbind, terms_dict), stringsAsFactors=FALSE)
colnames(terms2genes) <- "Genes"

#########################################################################
# Write output
#########################################################################

dir_out <- "data/human/GO/"
dir.create(dir_out, showWarnings=FALSE, recursive=TRUE)

outfn <- paste0(dir_out, "term_description.txt")
write.table(geo_df, outfn, row.names=TRUE, col.names=NA, sep="\t", quote=TRUE)

outfn <- paste0(dir_out, "genes2terms.txt")
write.table(genes2terms, outfn, row.names=TRUE, col.names=NA, sep="\t", quote=FALSE)

outfn <- paste0(dir_out, "terms2genes.txt")
write.table(terms2genes, outfn, row.names=TRUE, col.names=NA, sep="\t", quote=FALSE)

outfn <- paste0(dir_out, "gene_description.txt.txt")
write.table(genes_descr, outfn, row.names=TRUE, col.names=NA, sep="\t", quote=FALSE)

out.l <- list("Source" = "GO", 
              "Species" = "Human", 
              "genes2terms" = genes2terms, 
              "terms2genes" = terms2genes, 
              "term_description" = geo_df,
              "gene_description" = genes_descr)

outfn <- paste0(dir_out,  "pathway.rds")
saveRDS(out.l, outfn)

date_df <- date()
write.table(date_df, paste0(dir_out, "download_date.txt"), 
	row.names=FALSE, col.names=FALSE, quote=TRUE)

