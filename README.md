# snipets
code snipets with tricks and reusable components


#### Reading sparse matrix in python
```python
from scipy import io
counts = io.mmread('matrix.mtx.gz')
counts.toarray()
```
Tags : cellranger , count matrix, read sparse matrix

---

#### Getting HEX values from ggplot2 color hue - R code
```R
a <- scale_color_hue(h.start=180)
a$palette(7)#number of colors
```
This is helpful, especially if you are using both R and Python(or anythong else) and wanted to have plots have same color representation.

---

#### Exporting Seurat object from R to scanpy readable object

```R
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(reticulate)

libd.151507 <- readRDS("data/processed_data/151507.rds")#reading seurat object processed and stored in directory

SaveH5Seurat(libd.151507, filename = "data/processed_data/libd.151507.h5Seurat")
Convert("data/processed_data/libd.151507.h5Seurat", dest = "data/processed_data/libd.151507.h5ad")
```

```python
import pandas as pd
import numpy as np
import sys
from anndata import AnnData
import scanpy as sc

adata = sc.read_h5ad("data/processed_data/libd.151507.h5ad")

```
