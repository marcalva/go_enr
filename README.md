
# GO Enrichment

## Setup

## Downloading data

Run the following scripts to download ontology data from *GO*, *Reactome*, 
and *KEGG*.

```bash
Rscript R/get_go.R
Rscript R/get_kegg.R
Rscript R/get_go.human.R
```

Data is formatted into the following files

* **genes2terms.txt** text file that maps genes to ontology terms. The 
columns correspond to
    1. gene symbol
    2. comma-separated list of ontology term IDs.
* **terms2genes.txt** text file that maps ontology terms to genes in the 
term. The columns correspond to
    1. ontology term ID
    2. comma-separated list of gene symbols
* **term_description.txt** text file that gives the description of ontology 
terms. The first column gives the ontology term ID and later columns give 
descriptive information. There are always two columns: "Name" gives the 
short descripton and "Description" gives a more detailed description.
Additional meta data depends on the source of the annotations.
* **pathway.rds** RDS file that contains a list, where each of the 
elements gives a data frame. There are three data frames named 
"genes2terms", "terms2genes", and "term_description". There could be 
additional data frames depending on the annotation source.

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

You can also plot the results. The function `plot_enr` in the `R/plot.R` 
source file can be used to create a scatterplot of the odds ratio vs. the 
negative log 10 p-value of the enrichment terms.

```R
p <- plot_enr(go_df, color = "NameSpace")
```


