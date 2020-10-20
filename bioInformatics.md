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

Reading 10X spatial data in python

#Major part of follow code copied from 10X and other sources.

```python
import collections
import scipy.sparse as sp_sparse
import tables
 
CountMatrix = collections.namedtuple('CountMatrix', ['feature_ref', 'barcodes', 'matrix'])
 
def get_matrix_from_h5(filename):
    with tables.open_file(filename, 'r') as f:
        mat_group = f.get_node(f.root, 'matrix')
        barcodes = f.get_node(mat_group, 'barcodes').read()
        data = getattr(mat_group, 'data').read()
        indices = getattr(mat_group, 'indices').read()
        indptr = getattr(mat_group, 'indptr').read()
        shape = getattr(mat_group, 'shape').read()
        matrix = sp_sparse.csc_matrix((data, indices, indptr), shape=shape)
         
        feature_ref = {}
        feature_group = f.get_node(mat_group, 'features')
        feature_ids = getattr(feature_group, 'id').read()
        feature_names = getattr(feature_group, 'name').read()
        feature_types = getattr(feature_group, 'feature_type').read()
        feature_ref['id'] = feature_ids
        feature_ref['name'] = feature_names
        feature_ref['feature_type'] = feature_types
        tag_keys = getattr(feature_group, '_all_tag_keys').read()
        print("tag_keys:::",tag_keys)
        for key in tag_keys:
            feature_ref[key] = getattr(feature_group, 'genome').read()
         
        return CountMatrix(feature_ref, barcodes, matrix)
     
```
