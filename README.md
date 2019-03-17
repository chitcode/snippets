# snipets
code snipets with tricks and reusable components


#### Reading sparse matrix in python
```python
from scipy import io
counts = io.mmread('matrix.mtx.gz')
counts.toarray()
```
Tags : cellranger , count matrix

#### Getting HEX values from ggplot2 color hue - R code
```R
a <- scale_color_hue(direction = -1, h.start=180)
a$palette(7)#number of colors
```
