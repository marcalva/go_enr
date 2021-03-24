
# GO Enrichment

## Setup

These scripts can download ontology data from *GO*, *Reactome*, and 
*KEGG*. Data is downloaded and formatted into the following files

* **genes2terms.txt** text file where the first column gives the gene 
  name and the second gives a comma-separated list of IDs. These IDs 
  correspond to the annotation, ontology, or pathway of the resource.
* **terms2genes.txt** text file where the first column gives the 
  term ID and the second column gives a comma-separated list of gene 
  names. The terms IDs are the same as those listed above.
* **term_description.txt** text file where the first column gives 
  the term ID and later columns give meta data. 
  The column "Name" gives the short descrpition name of the ontology 
  term, while the column "Description" gives a longer description.
  Additional meta data depends on the source of the annotations.
* **pathway.rds** RDS file that contains each of the above data stored 
  in an RDS file. The RDS file contains a list of 3 elements, each 
  of the above.

## Downloading data

Run the following scripts

```bash
Rscript R/get_go.R
Rscript R/get_kegg.R
Rscript R/get_go.human.R
```

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


