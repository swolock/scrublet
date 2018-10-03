# Scrublet
**S**ingle-**C**ell **R**emover of Do**ublet**s  
  
Python code for identifying doublets in single-cell RNA-seq data. For details and validation of the method, see our preprint on [bioRxiv](https://www.biorxiv.org/content/early/2018/07/09/357368).

#### Quick start:
Given a raw counts matrix `counts_matrix` with cells as rows and genes as columns, calculate a doublet score for each cell with the following commands: 
```python
import scrublet as scr
scrub = scr.Scrublet(counts_matrix)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
```

Several optional parameters can influence the quality of the predictions. For a typical workflow, including interpretation of predicted doublet scores, see the example [notebook](./examples/scrublet_basics.ipynb).

#### Installation:
```bash
git clone https://github.com/AllonKleinLab/scrublet.git
cd scrublet
pip install -r requirements.txt
pip install --upgrade .
```
#### Other doublet detection tools:
[DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder)  
[DoubletDecon](https://github.com/EDePasquale/DoubletDecon)  
[DoubletDetection](https://github.com/JonathanShor/DoubletDetection)
