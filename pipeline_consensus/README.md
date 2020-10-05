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
### Linux
* create and activate the pipeline's environment
```bash
conda env create -f envs/covid.yml -n covid
conda activate covid
```
### MacOS
* install HTSlib from source using the [github repo here](https://github.com/samtools/htslib)
* install iVar from source using the [github repo here](https://github.com/andersen-lab/ivar)
* install GNU Utils for MacOS using Homebrew
```bash
brew install coreutils gnu-sed
```
* create and activate the pipeline's environment
```bash
conda env create -f envs/covid_macos.yml -n covid
conda activate covid
```

## Usage
* Ensure that the pipeline's Conda Environment is activated
* Specify user's data and parameters in `config.json` and save the file
* Run the pipeline using Snakemake
```bash
snakemake --cores [num_cores] --use-conda
```

## Development
In order to work on the development of the pipeline, the user must ensure that any additional dependencies are captured in the environment file `envs/covid.yaml`. Here are general instructions on how to do so:
* activate the pipeline's environment 
```bash
conda activate covid
```
* install any additional dependencies
```bash
conda install -c bioconda pokechop
```
* Update the environment YAML file with the additional dependencies
```bash
conda env export > envs/covid.yaml
```
* Update the remote repository with the new YAML file
```bash
git add envs/covid.yaml
git commit -m "added pokechop as a dependency"
```