# COVID Consensus Pipeline

## Pre-requisites
Prior to installing and running the pipeline, the following must be installed in the user's system

* Anaconda
```bash
wget https://repo.anaconda.com/archive/Anaconda3-2020.07-Linux-x86_64.sh
bash Anaconda3-2020.07-Linux-x86_64.sh
```

* iVar (version 1.2.2)

## Installation
* create and activate the pipeline's environment
```bash
conda env create -f envs/covid.yaml -n covid
conda activate covid
```

## Usage
* Ensure that the pipeline's Conda Environment is activated
* Specify user's data and parameters in `config.json` and save the file
* Run the pipeline using Snakemake
```bash
snakemake --cores [num_cores]
```