surface-based profiling of FCDs at 9.4T vs 3T

This is the code for the analyses in the manuscript
> Kronlage, ... Kuehn: "9.4 Tesla MRI in focal epilepsy patients with high-resolution surface-based profiling of focal cortical dysplasias", 2026 - in preparation.

# How to use:
## folder structure
- data in BIDS format is placed in `data/`, results and intermediate files will be stored in `data/derivatives`
- scripts are all in the root folder

## Requirements
Install conda environment 3T9T_profiling:
```
conda env create -f environment.yml
```
We need another environment for snakemake with a newer python version of 3.13
```
conda env create -f environment_snakemake.yml
```

Other dependencies:
- Freesurfer v 7.4.1. For reproducibility and to enable switching to version 7.2 (because of a bug in xhemi 7.4), the snakemake workflow can also use apptainer containers pulled from docker.
- laynii installed in ~/laynii/ (or change the path in plot_laynii.ipynb)

## Run freesurfer reconstructions and feature extraction using snakemake
```
./snakemake_local.sh 
```
... or run snakemake in the snakemake environment to specify the number of cores or other parameters, e.g., for parallelizing on a cluster

## Plots
### Figure 1-2
Notebook `plot_lesions_3T9T.ipynb``

### Figure 3
Script `plot_freeview_lesions_3T9T_surfs.py` for the MRI and surface overlays, notebook `plot_surf_features.ipynb` for the thickness features plots

### Figure 4
Notebook `plot_surf_features.ipynb`

### Figure 5
Notebook `plot_laynii.ipynb`


# Contact:
If you have any questions, comments or suggestions, get in touch:
<cornelius.kronlage@med.uni-tuebingen.de>