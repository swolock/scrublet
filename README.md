# Scrublet
**S**ingle-**C**ell **R**emover of Do**ublet**s  
  
Python code for identifying doublets in single-cell RNA-seq data. For details and validation of the method, see our preprint on [bioRxiv](https://www.biorxiv.org/content/early/2018/07/09/357368).

#### Quick start:
For a typical workflow, including interpretation of predicted doublet scores, see the example [notebook](./examples/scrublet_basics.ipynb).  
  
Given a raw (unnormalized) UMI counts matrix `counts_matrix` with cells as rows and genes as columns, calculate a doublet score for each cell: 
```python
import scrublet as scr
scrub = scr.Scrublet(counts_matrix)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
```
`scr.scrub_doublets()` simulates doublets from the observed data and uses a k-nearest-neighbor classifier to calculate a continuous `doublet_score` (between 0 and 1) for each transcriptome. The score is automatically thresholded to generate `predicted_doublets`, a boolean array that is `True` for predicted doublets and `False` otherwise. 

#### Installation:
```bash
git clone https://github.com/swolock/scrublet.git
cd scrublet
pip install -r requirements.txt
pip install --upgrade .
```
#### Other doublet detection tools:
[DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder)  
[DoubletDecon](https://github.com/EDePasquale/DoubletDecon)  
[DoubletDetection](https://github.com/JonathanShor/DoubletDetection)
