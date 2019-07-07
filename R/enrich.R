

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

# genes is a character vector of genes to test enrichment for
# background genes is a character vector of genes that are used as the background
# dir_go contains the GO files
# species is "human" or "mgi"
# min_genes only tests terms with at least this many genes
# method_correct corrects the FET p-value using this method. FDR by default
fet_go <- function(genes, 
				   GO_data, 
				   bg_genes=NULL, 
				   min_genes=30, 
				   method_correct="fdr"){
	# GO term annotations
	go_anno <- GO_data$go_anno
	# Gene annotations and their GO terms
	anno <- GO_data$anno
	# GO terms and their genes
	terms <- GO_data$terms

	# Set genes and background
	genes <- intersect(genes, rownames(anno))
	if (length(genes) == 0){
		stop("No genes given that are found in the GO annotation file")
	}
	if (is.null(bg_genes)){
		bg_genes <- rownames(anno)
	}
	bg_genes <- setdiff(bg_genes, genes)
	if (length(genes) == 0){
		stop("No background genes found in the GO annotation file")
	}

	# Loop through annotations
	enr <- list()
	for (term in rownames(terms)){
		term_genes <- terms[term,"Genes"]
		term_genes <- base::strsplit(x=term_genes, split=";")
		if (length(term_genes) < min_genes){
			next
		}
		# FET setup
		#  s1 | s2
		#  -------
		#  s3 | s4
		s1 <- base::intersect(genes, term_genes) # Genes in term
		s2 <- base::set_diff(genes, term_genes) # Genes NOT in term
		s3 <- base::intersect(bg_genes, term_genes) # BG in term
		s4 <- base::set_diff(bg_genes, term_genes) # BG NOT in term
		m <- matrix(c(length(s1), length(s2), length(s3), length(s4)), 
					nrow=2, ncol=2, byrow=TRUE)
		fet <- fisher.test(m)
		ret <- c(Name=go_anno[term, "Name"], 
				 NameSpace=go_anno[term, "NameSpace"], 
				 OR=fet$estimate, 
				 Null=fet$null.value, 
				 p=fet$p.value, 
				 GenesInTerm=length(s1), 
				 GenesOutTerm=length(s2), 
				 TermGenes=length(term_genes), 
				 Genes=s1)
		enr[[term]]
	}
	enr_df <- do.call(base::rbind, enr)
	enr_df$p_adj <- p.adjust(enr_df$p, method=method_correct)
	enr_df <- enr_df[order(enr_df$p_adj, decreasing=FALSE),]
	return(enr_df)
}

