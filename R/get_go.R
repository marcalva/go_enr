


dir_out = "data/ref/GO/"
dir.create(dir_out, showWarnings=FALSE, recursive=TRUE)



#########################################################################
# Gene ontology terms, definitions, descriptions.
#########################################################################

geo = readLines("http://snapshot.geneontology.org/ontology/go-basic.obo")
# geo is a list, each element is a line of the file go-basic.obo

geo_dict = list()
i=1
while (i <= length(geo)){
	if (geo[i] == "[Term]"){
		i = i+1
		id = sub("id: ", "", geo[i])
		i = i+1
		name = sub("name: ", "", geo[i])
		i = i+1
		namespace = sub("namespace: ", "", geo[i])
        i = i+1
        def = sub("def: ", "",  geo[i])
		geo_dict[[id]] = c(name, namespace, def)
	}
	i = i+1
}

geo_df = do.call(rbind, geo_dict)
colnames(geo_df) = c("Name", "NameSpace", "Definition")

write.table(geo_df, paste0(dir_out, "term_description.txt"),
			row.names=TRUE, col.names=NA, sep="\t", quote=FALSE)

#########################################################################
# Human
#########################################################################

#########################################################################
# Genes to terms and annotations
#########################################################################

con = gzcon(url("http://geneontology.org/gene-associations/goa_human.gaf.gz"))
txt = readLines(con)
goa = read.table(textConnection(txt), comment.char="!", header=FALSE, fill=TRUE, sep="\t",
				 stringsAsFactors=FALSE, quote="")

genes = unique(goa[,3])

goa_dict = list()
for (gene in genes){
	df = goa[goa[,3] == gene,]
	uniprots = paste(unique(df[,2]), collapse=",")
	descr = paste(unique(df[,10]), collapse=",")
	go_terms = paste(unique(df[,5]), collapse=",")
    goa_dict[[gene]] = c(go_terms, descr, uniprots)
}
goa_df = as.data.frame(do.call(rbind, goa_dict), stringsAsFactors=FALSE)
colnames(goa_df) = c("Terms", "Description", "UniProtID")

genes2terms <- goa_df[,"Terms",drop=FALSE]
outfn <- paste0(dir_out, "genes2terms.txt")
write.table(genes2terms, outfn, row.names=TRUE, col.names=NA, sep="\t", quote=FALSE)


colnames(goa_df) = c("UniProtID", "Description", "GO_Terms")

write.table(goa_df, paste0(dir_out, "GO_gene_terms.human.txt"),
			row.names=TRUE, col.names=NA, sep="\t", quote=FALSE)


#########################################################################
# Directly annotate
#########################################################################

mf = geo_df[geo_df[,"NameSpace"] == "molecular_function",]
bp = geo_df[geo_df[,"NameSpace"] == "biological_process",]
cc = geo_df[geo_df[,"NameSpace"] == "cellular_component",]

anno = as.data.frame(goa_df, stringsAsFactors=FALSE)

anno$MolecularFunction = sapply(anno$GO_Terms, function(x){
								terms = strsplit(x, ";")[[1]]
								terms = terms[terms %in% rownames(mf)]
								if (length(terms) == 0) NA
								else paste(mf[terms,"Name"], collapse=";")
			})

anno$BiologicalProcess = sapply(anno$GO_Terms, function(x){
								terms = strsplit(x, ";")[[1]]
								terms = terms[terms %in% rownames(bp)]
								if (length(terms) == 0) NA
								else paste(bp[terms,"Name"], collapse=";")
			})

anno$CellularComponent = sapply(anno$GO_Terms, function(x){
								terms = strsplit(x, ";")[[1]]
								terms = terms[terms %in% rownames(cc)]
								if (length(terms) == 0) NA
								else paste(cc[terms,"Name"], collapse=";")
			})

write.table(anno, paste0(dir_out, "GO_gene_terms.annotations.human.txt"),
			row.names=TRUE, col.names=NA, sep="\t", quote=FALSE)

#########################################################################
# Genes per ontology term
#########################################################################

terms = unique(goa[,5])
terms_dict = list()
for(term in terms){
	df = goa[goa[,5] == term,]
	genes = paste(unique(df[,3]), collapse=";")
	terms_dict[[term]] = genes
}

terms_df = as.data.frame(do.call(rbind, terms_dict), stringsAsFactors=FALSE)
colnames(terms_df) = "Genes"

write.table(terms_df, paste0(dir_out, "GO_terms.gene_list.human.txt"),
			row.names=TRUE, col.names=NA, sep="\t", quote=FALSE)









#########################################################################
# Mouse
#########################################################################
# Read in gene annotations
con = gzcon(url("http://geneontology.org/gene-associations/mgi.gaf.gz"))
txt = readLines(con)
goa = read.table(textConnection(txt), comment.char="!", header=FALSE, fill=TRUE, sep="\t",
				 stringsAsFactors=FALSE, quote="")

genes = unique(goa[,3])

goa_dict = list()
for (gene in genes){
	df = goa[goa[,3] == gene,]
	uniprots = paste(unique(df[,2]), collapse=";")
	descr = paste(unique(df[,10]), collapse=";")
	go_terms = paste(unique(df[,5]), collapse=";")
	goa_dict[[gene]] = c(uniprots, descr, go_terms)
}
goa_df = as.data.frame(do.call(rbind, goa_dict), stringsAsFactors=FALSE)
colnames(goa_df) = c("UniProtID", "Description", "GO_Terms")

write.table(goa_df, paste0(dir_out, "GO_gene_terms.mgi.txt"),
			row.names=TRUE, col.names=NA, sep="\t", quote=FALSE)

#########################################################################
# Directly annotate
#########################################################################

mf = geo_df[geo_df[,"NameSpace"] == "molecular_function",]
bp = geo_df[geo_df[,"NameSpace"] == "biological_process",]
cc = geo_df[geo_df[,"NameSpace"] == "cellular_component",]

anno = as.data.frame(goa_df, stringsAsFactors=FALSE)

anno$MolecularFunction = sapply(anno$GO_Terms, function(x){
								terms = strsplit(x, ";")[[1]]
								terms = terms[terms %in% rownames(mf)]
								if (length(terms) == 0) NA
								else paste(mf[terms,"Name"], collapse=";")
			})

anno$BiologicalProcess = sapply(anno$GO_Terms, function(x){
								terms = strsplit(x, ";")[[1]]
								terms = terms[terms %in% rownames(bp)]
								if (length(terms) == 0) NA
								else paste(bp[terms,"Name"], collapse=";")
			})

anno$CellularComponent = sapply(anno$GO_Terms, function(x){
								terms = strsplit(x, ";")[[1]]
								terms = terms[terms %in% rownames(cc)]
								if (length(terms) == 0) NA
								else paste(cc[terms,"Name"], collapse=";")
			})

write.table(anno, paste0(dir_out, "GO_gene_terms.annotations.mgi.txt"),
			row.names=TRUE, col.names=NA, sep="\t", quote=FALSE)

#########################################################################
# Genes per ontology term
#########################################################################

terms = unique(goa[,5])
terms_dict = list()
for(term in terms){
	df = goa[goa[,5] == term,]
	genes = paste(unique(df[,3]), collapse=";")
	terms_dict[[term]] = genes
}

terms_df = as.data.frame(do.call(rbind, terms_dict), stringsAsFactors=FALSE)
colnames(terms_df) = "Genes"

write.table(terms_df, paste0(dir_out, "GO_terms.gene_list.mgi.txt"),
			row.names=TRUE, col.names=NA, sep="\t", quote=FALSE)



# Save download date
date_df <- date()
write.table(date_df, paste0(dir_out, "download_date.txt"), 
	row.names=FALSE, col.names=FALSE, quote=TRUE)
