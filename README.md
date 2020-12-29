# Scrublet
**S**ingle-**C**ell **R**emover of Do**ublet**s  
  
Python code for identifying doublets in single-cell RNA-seq data. For details and validation of the method, see our paper in [Cell Systems](https://www.sciencedirect.com/science/article/pii/S2405471218304745) or the preprint on [bioRxiv](https://www.biorxiv.org/content/early/2018/07/09/357368).

#### Quick start:
For a typical workflow, including interpretation of predicted doublet scores, see the example [notebook](./examples/scrublet_basics.ipynb).  
  
Given a raw (unnormalized) UMI counts matrix `counts_matrix` with cells as rows and genes as columns, calculate a doublet score for each cell: 
```python
import scrublet as scr
scrub = scr.Scrublet(counts_matrix)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
```
`scr.scrub_doublets()` simulates doublets from the observed data and uses a k-nearest-neighbor classifier to calculate a continuous `doublet_score` (between 0 and 1) for each transcriptome. The score is automatically thresholded to generate `predicted_doublets`, a boolean array that is `True` for predicted doublets and `False` otherwise. 

#### Best practices:  
- When working with data from multiple samples, run Scrublet on each sample separately. Because Scrublet is designed to detect technical doublets formed by the random co-encapsulation of two cells, it may perform poorly on merged datasets where the cell type proportions are not representative of any single sample. 
- Check that the doublet score threshold is reasonable (in an ideal case, separating the two peaks of a bimodal simulated doublet score histogram, as in [this example](./examples/scrublet_basics.ipynb)), and adjust manually if necessary.
- Visualize the doublet predictions in a 2-D embedding (e.g., UMAP or t-SNE). Predicted doublets should mostly co-localize (possibly in multiple clusters). If they do not, you may need to adjust the doublet score threshold, or change the pre-processing parameters to better resolve the cell states present in your data.

#### Installation:
To install with PyPI:
```bash
pip install scrublet
```

To install from source:
```bash
git clone https://github.com/swolock/scrublet.git
cd scrublet
pip install -r requirements.txt
pip install --upgrade .
```

#### Old versions:
Previous versions can be found [here](./old_versions/).

#### Other doublet detection tools:
[DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder)  
[DoubletDecon](https://github.com/EDePasquale/DoubletDecon)  
[DoubletDetection](https://github.com/JonathanShor/DoubletDetection)
