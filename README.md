# snipets
code snipets with tricks and reusable components


#### Reading sparse matrix in python
```python
from scipy import io
counts = io.mmread('matrix.mtx.gz')
counts.toarray()
```
Tags : cellranger , count matrix
