# woublet
Python 2 code for identifying doublets in single-cell RNA-seq data.  
Alternatively, see https://github.com/swolock/scanpy for a [scanpy](https://github.com/theislab/scanpy) fork implementing this doublet detector. Example [here](./180601_scanpy_example.ipynb).

Given a raw counts matrix `E` in `scipy.sparse.csc_matrix` format, with cells as rows and genes as columns, you can calculate a doublet score for each cell with the following command: 
```python
import woublet
doublet_score = woublet.woublet(E)[0]
```

There are a number of optional parameters that can have major effects on the quality of results. For a typical workflow, including interpretation of predicted doublet scores, see the example [ipython notebook](./180306_basic_example.ipynb).

#### Libraries required to run woublet:
`numpy` (http://www.numpy.org/)  
`scipy` (https://www.scipy.org/)  
`sklearn` (http://scikit-learn.org/stable/)  

#### Additional recommended libraries:
`annoy` (https://github.com/spotify/annoy)  

#### Additional libraries required to run the example ipython notebook:
`matplotlib` (https://matplotlib.org/)  
`networkx` (https://networkx.github.io/)  
`fa2` (https://github.com/bhargavchippada/forceatlas2)  

#### Tested using these versions:
`python 2.7.12`  
`numpy 1.11.3`  
`scipy 0.19.0`  
`sklearn 0.18.1`  
`annoy 1.8.0`  
`matplotlib 1.5.1`  
`networkx 1.9.1`  
`fa2 0.2`  
