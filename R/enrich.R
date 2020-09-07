
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
		warning("No genes given that are found in the GO annotation file")
	}

	if (is.null(bg_genes)){
		bg_genes <- rownames(anno)
	} else { 
        bg_genes <- intersect(bg_genes, rownames(anno))
	    bg_genes <- setdiff(bg_genes, genes)
    }

	if (length(bg_genes) == 0){
		stop("No background genes found in the GO annotation file")
	}

	# Loop through annotations
	enr <- list()
	for (term in rownames(terms)){
		term_genes <- terms[term,"Genes"]
		term_genes <- base::strsplit(x=term_genes, split=";")[[1]]
		if (length(term_genes) < min_genes){
			next
		}
		# FET setup
		#  s1 | s2
		#  -------
		#  s3 | s4
		s1 <- base::intersect(genes, term_genes) # Genes in term
		s2 <- base::setdiff(genes, term_genes) # Genes NOT in term
		s3 <- base::intersect(bg_genes, term_genes) # BG in term
		s4 <- base::setdiff(bg_genes, term_genes) # BG NOT in term
		m <- matrix(c(length(s1), length(s2), length(s3), length(s4)), 
					nrow=2, ncol=2, byrow=TRUE)
		fet <- fisher.test(m)
		ret <- data.frame(Name=go_anno[term, "Name"], 
				 NameSpace=go_anno[term, "NameSpace"], 
				 OR=fet$estimate, 
				 Null=fet$null.value, 
				 p=fet$p.value, 
				 GenesInTerm=length(s1), 
				 GenesOutTerm=length(s2), 
				 TermGenes=length(term_genes), 
				 Genes=base::paste(s1, collapse=";"))
		enr[[term]] <- ret
	}
	enr_df <- do.call(base::rbind, enr)
	enr_df$p_adj <- p.adjust(enr_df$p, method=method_correct)
    enr_df$Term <- rownames(enr_df)
	enr_df <- enr_df[order(enr_df$p, decreasing=FALSE),]
	enr_df <- enr_df[,c("Term", "Name", "p", "p_adj",  "NameSpace", "OR", "Null", "GenesInTerm", "GenesOutTerm", "TermGenes", "Genes")]
	return(enr_df)
}

# genes is a character vector of genes to test enrichment for
# reactome is a list with genes2pathway and pathway2genes objects.
# bg_genes is a character vector of background genes (NULL by default).
# min_genes only tests terms with at least this many genes
# method_correct corrects the FET p-value using this method. FDR by default
fet_react <- function(genes, 
                      reactome, 
                      bg_genes=NULL, 
                      min_genes=30, 
                      method_correct="fdr"){
	# Genes to pathway
    react_genes <- rownames(reactome$genes2pathway)
    pathways <- rownames(reactome$pathway2genes)

	# Set genes and background
	genes <- intersect(genes, react_genes)
	if (length(genes) == 0){
		warning("None of genes are found in the annotations")
	}

	if (is.null(bg_genes)){
		bg_genes <- react_genes
	} else {
        bg_genes <- intersect(bg_genes, react_genes)
	    bg_genes <- setdiff(bg_genes, genes)
    }

	if (length(bg_genes) == 0){
		stop("No background genes found in the annotations")
	}

	# Loop through annotations
	enr <- list()
	for (pathway in pathways){
        path_genes <- reactome$pathway2genes[pathway,"Genes"]
		path_genes <- base::strsplit(x=path_genes, split=";")[[1]]
		if (length(path_genes) < min_genes){
			next
		}
		# FET setup
		#  s1 | s2
		#  -------
		#  s3 | s4
		s1 <- base::intersect(genes, path_genes) # Genes in pathway
		s2 <- base::setdiff(genes, path_genes) # Genes NOT in pathway
		s3 <- base::intersect(bg_genes, path_genes) # BG in pathway
		s4 <- base::setdiff(bg_genes, path_genes) # BG NOT in pathway
		m <- matrix(c(length(s1), length(s2), length(s3), length(s4)), 
					nrow=2, ncol=2, byrow=TRUE)
		fet <- fisher.test(m)
		ret <- data.frame(OR=fet$estimate, 
				          Null=fet$null.value, 
				          GenesInTerm=length(s1), 
				          GenesOutTerm=length(s2), 
				          TermGenes=length(path_genes), 
				          Genes=base::paste(s1, collapse=";"), 
				          p=fet$p.value)
		enr[[pathway]] <- ret
	}
	enr_df <- do.call(base::rbind, enr)
	enr_df$p_adj <- p.adjust(enr_df$p, method=method_correct)
    enr_df$Term <- rownames(enr_df)
	enr_df <- enr_df[order(enr_df$p, decreasing=FALSE),]
	enr_df <- enr_df[,c("Term", "OR", "Null", "GenesInTerm", "GenesOutTerm", "TermGenes", "Genes", "p", "p_adj")]
	return(enr_df)
}

