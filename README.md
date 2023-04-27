## Reproducibility repository for BENGAL ##

Author&Maintainer: Yuyao Song <ysong@ebi.ac.uk>

BENGAL: BENchmarking strateGies for cross-species integrAtion of singLe-cell RNA sequencing data

We provide the jupyter notebooks to reproduce the benchmarking results of this study. We also provide scripts to make the main figures. The starting materials would be the outputs by the BENGAL Nextflow pipeline. 

[Check out the BENGAL preprint](https://www.biorxiv.org/content/10.1101/2022.09.27.509674v1)
[Go to BENGAL pipeline main repo](https://github.com/Functional-Genomics/BENGAL)

### Contents:
1. Jupyter notebooks for benchmarking results ranking and aggregation are in notebooks/
2. Scripts for making figures are in scripts/
3. Final benchmarking metrics, scores (from different aggregation methods) and rankings are in results/

### How to run:
1. Create a conda environment with BENGAL_reproducibility.yml and activate this environment. Make sure the IRkernel is registered for jupyter sessions.
2. Open the jupyter notebook in a jupyter/juoyrer lab session, select the above kernel
3. Run cells in the notebook from top, all should run smoothly

#### License:
MIT license
