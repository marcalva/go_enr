
# Gene ontology enrichment

## Setup

The download scripts rely on the following packages

* [KEGGREST](https://bioconductor.org/packages/release/bioc/html/KEGGREST.html)
* [AnnotationDbi](https://bioconductor.org/packages/release/bioc/html/AnnotationDbi.html)
* [org.Hs.eg.db](http://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html)

The plotting functions require

* [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)
* [ggrepel](https://cran.r-project.org/web/packages/ggrepel/index.html)
* [scales](https://cran.r-project.org/web/packages/scales/index.html)
* [reshape2](https://cran.r-project.org/web/packages/reshape2/index.html)

## Downloading data

These scripts calculate enrichments for the gene ontology (**GO**), **Reactome**, 
and **KEGG** database annotations. The repo already comes with pre-downloaded 
data, but can be updated by downloading the latest versions with the 
following scripts.

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
terms. The first column gives the ontology term ID. The second column (Name) gives 
a descriptive name of the ontology term. The third column (Description) gives a 
desrciption of the term. The go file has a fourth column (NameSpace)
that classifies the ontology term into one of "biological_process", 
"cellular_component", or "molecular_function". The kegg file has a fourth 
column (Class) that gives the class of the ontology.
* **pathway.rds** RDS file that contains a list with three data frames 
corresponding to the files described above.
There could be additional data frames depending on the annotation source.

## Functions

The enrichment for pathways/ontologies is calculated using a Fisher's 
Exact Test (FET). The file `R/enrich.R` contains the enrichment functions
`fet_go` to perform this. Annotations can be read from an RDS object 
using the `readRDS` function.

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


