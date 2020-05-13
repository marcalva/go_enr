
# Read in GO data
read_GO <- function(dir_go="data/ref/GO/", species="human"){
	# GO term annotations
	go_anno_fn <- paste0(dir_go, "GO_id.name.namespace.txt")
	go_anno <- read.table(go_anno_fn, row.names=1, header=TRUE, stringsAsFactors=FALSE, sep="\t", comment="", quote="")

	# Read in annotations
	anno_fn <- paste0(dir_go, "GO_gene_terms.annotations.", species, ".txt")
	anno <- read.table(anno_fn, row.names=1, header=TRUE, stringsAsFactors=FALSE, sep="\t", comment="", quote="")

	# Read in GO terms and their genes
	terms_fn <- paste0(dir_go, "GO_terms.gene_list.", species, ".txt")
	terms <- read.table(terms_fn, row.names=1, header=TRUE, stringsAsFactors=FALSE, sep="\t", comment="", quote="")
		
	dd <- read.table(paste0(dir_go, "download_date.txt"), stringsAsFactors=FALSE)[1,1]

	ret <- list(go_anno=go_anno,
				anno=anno, 
				terms=terms, 
				date=dd)
	return(ret)
}

# Read in GO data
read_reactome <- function(dir_react="data/ref/Reactome/", species="human"){

	# Read in annotations
	g2p_fn <- paste0(dir_react, "genes2pathway.", species, ".txt")
	g2p <- read.table(g2p_fn, row.names=1, header=TRUE, quote="", 
                      stringsAsFactors=FALSE, sep="\t", comment="")

	# Read in GO terms and their genes
	p2g_fn <- paste0(dir_react, "pathway2genes.", species, ".txt")
	p2g <- read.table(p2g_fn, row.names=1, header=TRUE, stringsAsFactors=FALSE, sep="\t", comment="", quote="")
		
	dd <- read.table(paste0(dir_react, "download_date.txt"), 
                     stringsAsFactors=FALSE)[1,1]

	ret <- list(genes2pathway = g2p, 
                pathway2genes = p2g, 
				date=dd)
	return(ret)
}
