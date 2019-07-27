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
