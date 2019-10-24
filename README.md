
# GO Enrichment

## Downloading GO terms

Use the script:

```bash
Rscript R/get_go.R
```

to download the latest gene annotations from GO. They are 
downloaded to `data/ref/GO`. If you don't want to download them, you can 
use the ones already present downloaded July 6th 2019.

First, the gene ontologies themselves are downloaded to 
`data/ref/go/GO_id.name.namespace.txt`. This contains info on the 
ontology terms. The rownames are the GO ID (e.g. GO:0000001). The 
first column is the name of the GO term (e.g. mitochondrion inheritance) 
and the second column is the namespace (one of **biological_process**, 
**molecular_function** or **cellular component**). 

After downloading the ontologies, the genes are annotated with them. The 
annotations are downloaded to `data/ref/go/GO_gene_terms.human.txt` and 
`data/ref/go/GO_gene_terms.mgi.txt`. The row names are the gene IDs, the 
first column is the UniProtID, the second column is a description of the 
gene, and the third column contains the GO term IDs associated with gene. 
If there are multiple GO terms associated with the gene, they are separated 
by semicolons.

After getting the gene annotations, they are expanded so that we can directly 
read all the ontologies instead of the GO term ID. These expanded annotations 
are placed in the files `data/ref/go/GO_gene_terms.annotations.human.txt` and 
`data/ref/go/GO_gene_terms.annotations.mgi.txt`. The first columns are the 
same as the annotations: the row names are the gene IDs, the first column 
is the UniProtID, the second column is a description of the gene, and the 
third column contains the GO term IDs associated with gene. Then, the 
fourth, fifth, and sixth column contain the **MolecularFunction**, 
**BiologicalProcess**, and **CellularComponent** terms associated with the 
gene. Multiple terms are separated by semicolons.

Then, the reverse is done where for each GO term ID, all the gene names associated 
with it are listed in the first column. The files are in 
`data/ref/go/GO_terms.gene_list.human.txt` and 
`data/ref/go/GO_terms.gene_list.mgi.txt`.

## Functions

With the GO data downloaded above, you can use the functions to test for 
enrichment based on a gene set. To read in the anntoations, use the 
function `read_GO` in the `R/read.R` file.

The function `fet_go` is used to test enrichment of a set of genes using 
a Fisher's Exact Test (FET). This function is in the `R/enrich.R`.

```R
go_dat <- read_GO(dir_go="data/ref/GO/", species="human")
go_df <- fet_go(gene_set, bg_genes=NULL, GO_data=go_dat)
```
