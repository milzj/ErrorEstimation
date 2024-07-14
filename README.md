[![GNU GPLv3 License](https://img.shields.io/github/license/milzj/FW4PDE)](https://choosealicense.com/licenses/gpl-3.0/)
[![Test](https://github.com/milzj/FW4PDE/actions/workflows/test-FW4PDE.yml/badge.svg?style=plastic)](https://github.com/milzj/ErrorEstimation/actions/workflows/test.yml)
[![arXiv](https://img.shields.io/badge/arXiv-2402.15948-b31b1b.svg)](http://arxiv.org/abs/2402.15948)
[![DOI](https://zenodo.org/badge/751907933.svg)](https://zenodo.org/doi/10.5281/zenodo.10798891)


# Supplementary code for the manuscript: Criticality Measure-based Error Estimates for Infinite Dimensional Optimization

This repository contains supplementary code for the manuscript
> Danlin Li, Johannes Milz 2024
> Criticality Measure-based Error Estimates For Infinite Dimensional Optimization, preprint, http://arxiv.org/abs/2402.15948

## Abstract

Motivated by optimization with differential equations, we consider optimization problems with Hilbert spaces as decision spaces. As a consequence of their infinite dimensionality, the numerical solution necessitates finite dimensional approximations and discretizations. We develop an approximation framework and demonstrate criticality measure-based error estimates. We consider criticality measures inspired by those used within optimization methods, such as semismooth Newton and (conditional) gradient methods. Furthermore, we show that our error estimates are order-optimal. Our findings augment existing distance-based error estimates, but do not rely on strong convexity or second-order sufficient optimality conditions. Moreover, our error estimates can be used for code verification and validation. We illustrate our theoretical convergence rates on linear, semilinear and bilinear PDE-constrained optimization.


## Reproducing the numerical simulations
All the scripts are located in the folder called `code` in the repository. Is is assumed that you run the script from within this folder.

### Installation

```
conda env create -f environment.yml
conda activate ErrorEstimation
```
### Simulation

We can solve the optimization problems and evaluate criticality measures by running the shell script.

```
./simulation_rates.sh
```
This will create a new folder in `code/output`. Its name is a time stamp.

### Plotting

[After updating the time stamp](https://github.com/milzj/ErrorEstimation/blob/79611618e38d88e5627e6b37275bf6d2c62bbfe3/code/plot_rates.sh#L1), the convergence rate plots can be created by running

```
./plot_rates.sh
```

[After updating the time stamp](https://github.com/milzj/ErrorEstimation/blob/79611618e38d88e5627e6b37275bf6d2c62bbfe3/code/plot_control.sh#L1), 
the optimal controls can be visualized by running

```
./plot_control.sh
```

## Using Docker image

The Dockerfile can be used by running:

```
docker build -t errorestimation . --no-cache --network=host
docker run -it errorestimation
```

The prebuild docker image can be downloaded and run by executing

```
docker pull ghcr.io/milzj/errorestimation:latest
docker run --rm -it ghcr.io/milzj/errorestimation:latest
```
The latter command my be replaced by

```
docker run -ti -v ${PWD}:/root/shared -w /root/shared --entrypoint=/bin/bash --rm ghcr.io/milzj/errorestimation:latest
```

## Citation

```
@software{Li2024,
  author = {Li, Danlin and Milz, Johannes},
  doi = {10.48550/arXiv.2402.15948},
  title = {Supplementary code for the manuscript: Criticality Measure-based Error Estimates for Infinite Dimensional Optimization},
  url = {https://github.com/milzj/ErrorEstimation},
  year = {2024}
}
```

## Having issues
If you have any troubles please file an issue in the GitHub repository.

## Authors

- Danlin Li
- Johannes Milz 
