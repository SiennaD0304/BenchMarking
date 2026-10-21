# BenchMarking 
Orthologs are a type of homolog which occur due to speciation, resulting in conserved sequences and function. As a result, orthologs are frequently used in comparative genomics, and gene annotation, and determining species phylogeny. The use of orthologs in bioinformatics has resulted in the development of numerous tools. In particular, recent influx of whole genome sequencing has resulted in a greater need for fast, and accurate ortholog tools. Here we present benchmarking of four recently made, or updated tools: OrthoFinder V3, Broccoli V2, FastOMA and OrthoPhy. 


# Installation
To ensure that our results are replicable, instructions on how to install each tool are given below. 

First start by installing Mamba if not already installed to help increase the speed of making environments in the next steps.
```
# installing mamba 
conda install -n base -c conda-forge mamba
```
If Mamba is not supported or unable to be installed all environments can be created using Conda.

Once Mamba is installed clone the repo below. 
```
git clone https://github.com/SiennaD0304/BenchMarking.git 
cd ortho-benchmark
```

Create Conda environments:
The environments are needed to run each of the ortholog tools due to the large amount of dependancies associated with each tool. These environments are all created using .yml files. The trees environment allows for the visualization of the gene trees outputted by each tool. 
```
# creating conda environments
mamba env create -f environments/orthofinder_environment.yml
mamba env create -f environments/fastoma_environment.yml
mamba env create -f environments/broccoli_environment.yml
mamba env create -f environments/orthophy_environment.yml
mamba env create -f environments/trees_environment.yml
```

Installing databases and dependancies:
OrthoFinder, Broccoli and FastOMA are all available through bioconda. However Orthophy is not and therefore their Github must be installed directly. 
```
mkdir -p tools
git clone https://github.com/tmakwtnb/OrthoPhy.git tools/OrthoPhy 
```
We will then confirm that all of our environments were successfully installed using the command below. 
```
conda env list
```

The tool FastOMA requires the OMA Orthology Database to be installed to successfully run. This database must be installed manually as the file size is too large to be stored in Github. To ensure the database is accessible, keep in the data directory. Note: This is a large database and will therefore take about an hour to fully download.

```
cd data
wget https://omabrowser.org/All/LUCA.h5 &
```

# Testing Tools with Practice Dataset

To verify if the tools were installed 
```
conda run -n of3 orthofinder -h | head -5
conda run -n broccoli broccoli --help | head -5
conda run -n fastoma nextflow -version
conda run -n fastoma omamer --help | head -3
conda run -n orthophy python tools/OrthoPhy/OrthoPhy.py --help | head -10
conda run -n treepipe iqtree2 --version
```

To verify if the tools run on test data and ensure that all dependancies were installed correctly, the python script test_tools.py was used. This script checks that all directories present contain a file, and produce an output using a series of pass fail tests on a very small set of known protein orthologs. You must be be in the main folder to run the following commands 
```
#email is needed because fasta files are being pulled from NCBI
python test/test-scripts/get_test_data.py -e user@email.com

python test/test-scripts/prep_inputs.py -i test/test-data/test_raw -o test/test-data/prepped

#this command will test all of the tools
python test/test-scripts/test_tools.py -i test/test-data/prepped -o results/test-data 

#if the previous command fails you can test individual tools as well
python test/test-scripts/test_tools.py -i test/test-data/prepped -o results/test-data -t (insert tool name here)
```

# Results of Testing Data
A successful run of the script test_tools.py will allow for OrthoFinder, Broccoli and FastOMA to pass all tests. However, OrthoPhy will not pass. This is due to the very small sample size used for the test. OrthoPhy normalizes ortholog groups to increase the tool's speed, causing the very small sample size to result in no output. However, because OrthoPhy runs as an updated version of OrthoFinder if OrthoFinder passes then OrthoPhy will as well with a larger dataset. 





