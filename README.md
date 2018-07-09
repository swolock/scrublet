# Scrublet
**S**ingle-**C**ell **R**emover of Do**ublet**s  
Python code for identifying doublets in single-cell RNA-seq data.  
Alternatively, see https://github.com/swolock/scanpy for a [scanpy](https://github.com/theislab/scanpy) fork implementing this doublet detector. Example [here](./180601_scanpy_example.ipynb).

Given a raw counts matrix `E` in `scipy.sparse.csc_matrix` format, with cells as rows and genes as columns, you can calculate a doublet score for each cell with the following command: 
```python
import scrublet as scr
doublet_score = scr.predict_doublets(E)[0]
```

There are a number of optional parameters that can have major effects on the quality of results. For a typical workflow, including interpretation of predicted doublet scores, see the example [ipython notebook](./180306_basic_example.ipynb).

#### To install:
```bash
git clone https://github.com/AllonKleinLab/scrublet.git
cd scrublet
pip install -r requirements.txt
pip install --upgrade .
```
