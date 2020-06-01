## Converting Gene symbols to gene id (ENTREZID or NCBI etc) format

```r
library(clusterProfiler)
data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]
gene.df <- bitr(gene, fromType = "ENTREZID",
        toType = c("ENSEMBL", "SYMBOL"),
        OrgDb = org.Hs.eg.db)
head(gene.df)
```
source : https://yulab-smu.github.io/clusterProfiler-book/chapter5.html

There are other alternative approches for these conversions.

-----

Getting genes list for GO term

Ref : https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.html
```r
library("biomaRt")
#Mart for human genes
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

genes.df <- getBM(attributes=c('hgnc_symbol', 'go_id'),
+       filters = 'go', values = 'GO:0008009', mart = ensembl)

# Somewhat this is not filtering based on the GO value
go.genes <- genes.df$hgnc_symbol[df$go_id == "GO:0008009"]
```
----
