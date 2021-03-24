
#' Enrichment test using Fisher's Exact test
#' 
#' Takes genes in @p genes and test for enrichment in the terms listed in 
#' @p terms2genes. The 2-by-2 contingency table is formed by taking the 
#' query genes in the pathway, the query genes not in the pathway, the 
#' background genes in the pathway, and the pathway genes not in the pathway. 
#' The query genes are remoed from the background genes before forming 
#' the contingency table.
#' 
#' @param genes character vector of query genes of interest to test 
#'  enrichment for
#' @param genes2terms data frame where the row names correspond to genes.
#' @param terms2genes data frame where the row names correspond to ontology terms 
#'  and the first column contains a comma-separated list of gene symbols contained 
#'  in the ontology term.
#' @param bg_genes list of background genes for the NULL.
#' @param min_genes minimum number of genes in an ontology term to test 
#'  enrichment for.
#' @param prior integer count, which is added to all counts in the 2-by-2 
#'  contingency table before running Fisher's Exact test.
#' @param method_correct p-value correction method
#' @return data frame with ontology enrichments, sorted by p-value.
fet_enr <- function(genes, 
                    genes2terms, 
                    terms2genes, 
                    bg_genes = NULL, 
                    min_genes = 30, 
                    prior = 1, 
                    method_correct = "fdr"){
	# Genes to pathway
    all_genes <- rownames(genes2terms)
    all_terms <- rownames(terms2genes)

    terms2genes[1,] <- as.character(terms2genes[1,])

	# Set genes and background
	genes <- intersect(genes, all_genes)
	if (length(genes) == 0){
		warning("None of genes are found in the annotations")
	}

	if (is.null(bg_genes)){
		bg_genes <- all_genes
	}

	if (length(bg_genes) == 0){
		stop("No background genes found in the annotations")
	}

    bg_genes <- setdiff(bg_genes, genes)
    # print(length(bg_genes))
    if (length(bg_genes) < 30L){
        warning("only ", as.character(length(bg_genes)), " background genes")
    }

	# Loop through annotations
	enr <- list()
	for (p in all_terms){
        p_genes <- terms2genes[p,1]
        p_genes <- strsplit(x = p_genes, split = ",")[[1]]
		if (length(p_genes) < min_genes){
			next
		}
		# FET setup
		#  s1 | s2
		#  -------
		#  s3 | s4
		s1 <- base::intersect(genes, p_genes) # Genes in pathway
		s2 <- base::setdiff(genes, p_genes) # Genes NOT in pathway
		s3 <- base::intersect(bg_genes, p_genes) # BG in pathway
		s4 <- base::setdiff(bg_genes, p_genes) # BG NOT in pathway
		m <- matrix(c(length(s1) + prior, 
                      length(s2) + prior, 
                      length(s3) + prior,
                      length(s4) + prior), 
					nrow=2, ncol=2, byrow=TRUE)
		fet <- fisher.test(m)
		ret <- data.frame(OR=fet$estimate, 
				          Null=fet$null.value, 
				          GenesInTerm=length(s1), 
				          GenesOutTerm=length(s2), 
				          TermGenes=length(p_genes), 
				          Genes=base::paste(s1, collapse=","), 
				          p=fet$p.value)
		enr[[p]] <- ret
	}
	enr_df <- as.data.frame(do.call(base::rbind, enr))
	enr_df$p_adj <- p.adjust(enr_df$p, method=method_correct)
	enr_df <- enr_df[order(enr_df$p, decreasing=FALSE),]
	enr_df <- enr_df[,c("OR", "Null", "GenesInTerm", "GenesOutTerm", "TermGenes", "Genes", "p", "p_adj")]
	return(enr_df)
}

