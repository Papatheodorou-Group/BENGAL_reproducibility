## Reproducibility repository for BENGAL ##

Author&Maintainer: Yuyao Song <ysong@ebi.ac.uk>

BENGAL: BENchmarking strateGies for cross-species integrAtion of singLe-cell RNA sequencing data

We provide the jupyter notebooks and scripts to reproduce the benchmarking results of the BENGAL study. 

The starting materials would be the outputs by the BENGAL Nextflow pipeline. The same notebooks could be used for other benchmarking runs, while here we also provide combined results in .tsv format from the study.

[Check out the BENGAL preprint](https://www.biorxiv.org/content/10.1101/2022.09.27.509674)

[Go to BENGAL pipeline main repo](https://github.com/Functional-Genomics/BENGAL)

### Contents:
1. Jupyter notebooks for benchmarking results ranking and aggregation are in notebooks/
2. Scripts for making figures are in scripts/
3. Final benchmarking metrics, scores and rankings are in results/

### How to run:
1. Create a conda environment with BENGAL_reproducibility.yml. 
2. Register the IRkernel in the above environment into jupyter sessions. Refer to [IRkernel](https://github.com/IRkernel/IRkernel).
2. Open the jupyter notebook in a jupyter/jupyter lab session, select the above kernel.
3. Run cells in the notebook from top, all should run smoothly

#### License:
MIT license
