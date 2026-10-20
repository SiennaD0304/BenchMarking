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
git clone https://github.com/SiennaD0304/BenchMarking.git 
cd ortho-benchmark
```
Create Conda environments. The environments are needed to run each of the ortholog tools. These environments are all created using .yml files. The trees environment allows for the visualization of the gene trees outputted by each tool. 
```
# creating conda environments
mamba env create -f z/environments/orthofinder_environment.yml
mamba env create -f z/environments/fastoma_environment.yml
mamba env create -f z/environments/broccoli_environment.yml
mamba env create -f z/environments/orthophy_environment.yml
mamba env create -f z/environments/trees_environment.yml
```

OrthoFinder, Broccoli and FastOMA are all available through bioconda. However Orthophy is not and therefore their Github must be installed directly. 
```
mkdir -p tools git clone https://github.com/tmakwtnb/OrthoPhy.git tools/OrthoPhy 
```
We will then confirm that all of our environments were successfully installed using the command below. 
```
conda env list
```
# Installing databases and dependancies 
The tools FastOMA required a specific database to be installed to successfully run. To store their database make you are in the data folder available through our Github. Note: This is a large database and will therefore take about an hour to fully download
```
cd data
wget https://omabrowser.org/All/LUCA.h5 &
```

# Testing Tools 
```
# verify each tool runs
conda run -n of3 orthofinder -h | head -5
conda run -n broccoli broccoli --help | head -5
conda run -n fastoma nextflow -version
conda run -n fastoma omamer --help | head -3
conda run -n orthophy python tools/OrthoPhy/OrthoPhy.py --help | head -10
conda run -n treepipe iqtree2 --version
``` 
