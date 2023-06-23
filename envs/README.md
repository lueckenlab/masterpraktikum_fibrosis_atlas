# Conda environments

## Getting started

Install some base environments for using scanpy and/or R packages.

Use [mamba](https://mamba.readthedocs.io/en/latest/) (faster) or conda to create a new environment from a YAML file.

```shell
mamba env create -f envs/scanpy.yaml
```

### Get started with Jupyter Lab

```shell
mamba env create -f envs/jupyterlab.yaml
conda activate jupyterlab
jupyter lab
```

As long as environments have `ipykernel` and/or `r-irkernel` installed, you can access the environments from the jupyter lab interface.

### Install other environments

You can install the other environments in order to use them:

```shell
mamba env create -f envs/scanpy.yaml
mamba env create -f envs/scanpy_r.yaml
```

### Update an environment

If you want to update an environment that you have already installed, you can use the `update` command.
This applies if you want to add more dependencies to your environment.
Adding it to the environment file ensures that your workflow stays reproducible and dependencies get documented accordingly.

```shell
mamba env update -f envs/scanpy.yaml
```

## Create your own environment

If you want to use non-standard tools or specific workflows that can be hard to install, we recommend you to create a new environment.
You can use the existing environments as a template.

e.g.
```yaml
name: my_env_name
channels:
  - conda-forge
  - bioconda
dependencies:
  - python=3.9
  - scanpy=1.9
  - numpy<1.24
  - matplotlib<3.7
  - ipykernel
  - leidenalg
  # add more dependencies
```

Give your environment a meaningful name. Having the same name for the file and the environment makes things easier.

```shell
mamba env create -f envs/my_env_name.yaml
```



