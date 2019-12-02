
library(AnnotationDbi)
library(org.Hs.eg.db)


dir_out = "data/ref/Reactome/"
dir.create(dir_out, showWarnings=FALSE, recursive=TRUE)

# Reactome database
rea <- readLines('https://reactome.org/download/current/NCBI2Reactome_All_Levels.txt')
rea_l <- lapply(rea, function(i) strsplit(i, "\t")[[1]])
rea_df <- do.call(rbind, rea_l)
rea_df <- rea_df[rea_df[,6] == "Homo sapiens",]

# symbols <- mapIds(org.Hs.eg.db, rea_df[,1], 'SYMBOL', 'ENSEMBL')
symbols <- mapIds(org.Hs.eg.db, rea_df[,1], 'SYMBOL', 'ENTREZID')
rea_df <- cbind(symbols, rea_df)
rea_df <- rea_df[!is.na(rea_df[,1]),]

genes <- unique(rea_df[,1])

pathways <- unique(rea_df[,5])

rea_dict <- list()
pathways_dict <- list()
for (gene in genes){
    datf <- rea_df[rea_df[,1] == gene,,drop=FALSE]
    hsa <- datf[,3]
    p_url <- datf[,4]
    descr <- datf[,5]
    keep <- !duplicated(descr)
    hsa <- paste(hsa[keep], collapse=";")
    p_url <- paste(p_url[keep], collapse=";")
    descr <- paste(descr[keep], collapse=";")
    rea_dict[[gene]] <- c(hsa, 
                          p_url, 
                          descr)
    for (p in datf[,5]){
        pathways_dict[[p]] <- c(pathways_dict[[p]], gene)
    }
}

for (p in names(pathways_dict)){
    g <- unique(pathways_dict[[p]])
    g <- paste(g, collapse=";")
    pathways_dict[[p]] <- g
}

rea_datf <- as.data.frame(do.call(rbind, rea_dict), stringsAsFactors=FALSE)
write.table(rea_datf, 
            paste0(dir_out, "genes2pathway.human.txt"),
			row.names=TRUE, 
            col.names=NA, 
            sep="\t", 
            quote=FALSE)


path_datf <- as.data.frame(do.call(rbind, pathways_dict), stringsAsFactors=FALSE)
colnames(path_datf) <- "Genes"
write.table(path_datf, 
            paste0(dir_out, "pathway2genes.human.txt"),
			row.names=TRUE, 
            col.names=NA, 
            sep="\t", 
            quote=FALSE)

reactome <- list("genes2pathway" = rea_datf, 
                 "pathway2genes" = path_datf)
saveRDS(reactome, paste0(dir_out, "Reactome.RDS"))

# Save download date
date_df <- date()
write.table(date_df, 
            paste0(dir_out, "download_date.txt"), 
            row.names=FALSE, 
            col.names=FALSE, 
            quote=TRUE)

