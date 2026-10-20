# BenchMarking 
Orthologs are a type of homolog which occur due to speciation, resulting in conserved sequences and function. As a result, orthologs are frequently used in comparative genomics, and gene annotation, and determining species phylogeny. The use of orthologs in bioinformatics has resulted in the development of numerous tools. In particular, recent influx of whole genome sequencing has resulted in a greater need for fast, and accurate ortholog tools. Here we present benchmarking of four recently made, or updated tools: OrthoFinder V3, Broccoli V2, FastOMA and OrthoPhy. 


# Installation
to ensure that our results are replicable, instructions on how to install each tool are given below. 

First start by installing mamba if not already installed to help increase the speed of making environments in the next steps.
```
# installing mamba 
conda install -n base -c conda-forge mamba
```
If mamba is not supported or unable to be installed all environments can be created only using Conda.


Clone Repo: 
```
git clone [REPO_URL]
cd ortho-benchmark
```
Create Conda environments. The environments are needed to run each of the ortholog tools. These environments are all created using .yml files. 
```
# creating conda environments
mamba env create -f z/environments/orthofinder_environment.yml
mamba env create -f z/environments/fastoma_environment.yml
mamba env create -f z/environments/broccoli_environment.yml
mamba env create -f z/environments/orthophy_environment.yml
mamba env create -f z/environments/trees_environment.yml
```


