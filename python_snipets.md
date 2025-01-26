#### Reading sparse matrix in python
Helpful for reading scRNAseq data.
```python
from scipy import io
counts = io.mmread('matrix.mtx.gz')
```
Tags : cellranger , count matrix

---
#### One-hot-encoding in pure python

```python
def convert_one_hot_encode(y,classes):
    Y = np.eye(classes)[y.reshape(-1).T]
    return Y
```
---
#### Install conda package from Jupyter notebook cell
Prob: Conda needs permission before installing any package, however when you are installing any package from Jupyter notebook cell, you cannot give permission interactively.
```python
import sys
!conda install --yes --prefix {sys.prefix} numpy
```

---

#### Conda new environment

Prob: After creating a new environment in conda, to fix Jupyter notebooks installation issues.
```sh
conda create --name newenv
conda activate newenv
conda install nb_conda_kernels
```

#### Google colab disconnects
Google colab disconnects if the browser is not active for a while.
To avoid disconnection, execute the following code in the chrome console. Inspired from the [blogpost](https://medium.com/@shivamrawat_756/how-to-prevent-google-colab-from-disconnecting-717b88a128c0).
It seems the code in the blog has broken. So I modified to click on the files button in regular interval.
```javascript
function ClickConnect(){
console.log("Working");
files_btn = document.querySelector("body > div.notebook-vertical.colab-left-pane-open > div.notebook-horizontal > colab-left-pane > div > paper-listbox > paper-item:nth-child(3) > paper-icon-button").shadowRoot.querySelector("#icon");
files_btn.click()
}
setInterval(ClickConnect,60000)
```

#### Increase matplotlib image quality (dpi) in Jupyter notebook
```python
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
```
