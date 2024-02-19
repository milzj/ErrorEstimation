[![GNU GPLv3 License](https://img.shields.io/github/license/milzj/FW4PDE)](https://choosealicense.com/licenses/gpl-3.0/)
[![Test](https://github.com/milzj/FW4PDE/actions/workflows/test-FW4PDE.yml/badge.svg?style=plastic)](https://github.com/milzj/ErrorEstimation/actions/workflows/test.yml)
[![arXiv](https://img.shields.io/badge/arXiv-2306.17032-b31b1b.svg)](https://arxiv.org/abs/2306.17032)

# Supplementary code for the manuscript: Criticality Measure-based Error Estimates For Infinite Dimensional Optimization

This repository contains supplementary code for the manuscript
> A. Author, B. Author Year
> Criticality Measure-based Error Estimates For Infinite Dimensional Optimization, Journal of ..., volume, page, url

## Abstract

## Installation

```
conda env create -f environment.yml
conda activate ErrorEstimation
```

## Reproducing the numerical simulations
All the scripts are located in the folder called `code` in the repository. Is is assumed that you run the script from within this folder.

### Pre-processing
In order to reproduce the results you need to first run the pre-processing script
```
python3 pre_processing.py
```
This will convert the meshes from Gmsh to a dolfin format.

### Simulation
The next step is to run the fiber generation. You can do this by running the script
```
python3 run_fiber_generation.py
```
This will create a new folder `code/results` containing files called `microstructure_<heart_nr>.h5`.

### Plotting
The final step is to postprocess the results by running the script
```
python3 postprocess.py
```
This will generate a file for visualizing the fibers in the Paraview (inside `code/results` called  `fiber_<heart_nr>.xdmf`). This script will also compare some features computed from the fibers with the results published in the (artificial) paper. If the results differ, then the program will raise an error.

## Citation

```
@software{Lisa_My_Research_Software_2017,
  author = {Lisa, Mona and Bot, Hew},
  doi = {10.5281/zenodo.1234},
  month = {12},
  title = {{My Research Software}},
  url = {https://github.com/scientificcomputing/example-paper},
  version = {2.0.4},
  year = {2017}
}
```

## Having issues
If you have any troubles please file and issue in the GitHub repository.

## Authors

- A. Author
- B. Author
