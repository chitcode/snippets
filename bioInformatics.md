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

--

## Getting genes list for a GO term

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
--

## Reading 10X spatial data in python

*Major part of the following code taken from 10X and other sources.*

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
--

##Reading 10X Visium Spatial data 

For reading 10X Visium outputs, Load10X\_Spatial function from Seurat package is very convenient. Often times we download data from published papers, but the data repositories do not contain the filtered\_feature\_bc\_matrix.h5 file, which is a required file for Load10X_Spatial. The below code will help to generat the filtered\_feature\_bc\_matrix.h5 file so that Load10X\_Spatial function can be used.

Note that, after downloading the data from the repositories, you need to place them in the required folder structure.

```
path/to/data/dir
|__filtered_feature_bc_matrix
|   |__ barcodes.tsv.gz
|   |__ features.tsv.gz
|   |__ matrix.mtx.gz
|__spatial
    |__ tissue_lowres_image.png
    |__ scalefactors_json.json
    |__ tissue_positions_list.csv
    
```


```r
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#
#BiocManager::install("DropletUtils")

library(DropletUtils)
library(Seurat)

write10X_to_h5 <- function(path){
  coount.data <- readMM(paste0(path,"filtered_feature_bc_matrix/matrix.mtx.gz"))
  barcodes <- read.table(paste0(path,"filtered_feature_bc_matrix/barcodes.tsv.gz"),col.names="barcode")
  
  features <- read.table(paste0(path,"filtered_feature_bc_matrix/features.tsv.gz"), col.names = c("id","name","feature_type","expr"))
  write10xCounts(paste0(path,"/filtered_feature_bc_matrix.h5"), coount.data, gene.id=features$id, 
    gene.symbol=features$name, barcodes=barcodes$barcode, type="HDF5")
}

write10X_to_h5("path/to/data/dir")

spatial_data <- Load10X_Spatial("path/to/data/dir")
```
--

## Getting genes for a KEGG pathways

Here we are looking for genes in Melanoma [Pathway](https://www.genome.jp/kegg-bin/get_htext?query=05218&htext=br08901.keg&option=-a&node_proc=br08901_org&proc_enabled=hsa).

```r
library(KEGGREST)
keggReturns <- keggGet("hsa05218")[[1]]$GENE

#genes are at alternate lines, so first we will extract the lines containing gene symbols
genesLines <- seq(from=2, to = length(keggReturns), by = 2)
keggReturns <- keggReturns[genesLines]
#further cleaning
genesList <- gsub("\\;.*","",keggReturns)
head(genesList)
#"FGF1"  "FGF2"  "FGF3"  "FGF4"  "FGF17" "FGF6" 
```
