# Conda environments

## Getting started

Install some base environments for using scanpy and/or R packages.

Use [mamba](https://mamba.readthedocs.io/en/latest/) (faster) or conda to create a new environment from a YAML file.

```
mamba env create -f envs/scanpy.yaml
```

Get started with Jupyter Lab

```
mamba env create -f envs/jupyterlab.yaml
conda activate jupyterlab
jupyter lab
```

As long as environments have `ipykernel` and/or `r-irkernel` installed, you can access the environments from the jupyter lab interface.

