#### Reading sparse matrix in python

```python
from scipy import io
counts = io.mmread('matrix.mtx.gz')
```
Tags : cellranger , count matrix


---
#### Install conda package from Jupyter notebook cell
```python
import sys
!conda install --yes --prefix {sys.prefix} numpy
```

---

#### Conda new environment
```sh
conda create --new newenv
conda install nb_conda_kernels
```
