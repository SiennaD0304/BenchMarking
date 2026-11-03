#Species trees are required for the Tool FastOMA to run. This tree must be in Newick formatting and consist 
# of each species present in the testing dataset. The Quest for Ortholog database provides a list of taxaIDs 
# for all species present in their test dataset from their README, this script translates the IDs to species and forms a species
#  tree in Newick format 

#to run this script using argparse, run the following commands in terminal 
# cd BenchMarking/TestData/Testing_datasets
#python3 making_newick_tree.py --input TaxonomyIDs_2020 --output tree.txt

import argparse 
import sys
import numpy as np

#used to parse command line argument 
def check_arg(args=None):
   parser = argparse.ArgumentParser(description="Produce a newick species tree for FastOMA")
   parser.add_argument("-i", "--input",
         help = "path to input file", 
         required = "True"
         )
   parser.add_argument("-o", "--output", 
         help = "output file name", 
         required = "True"
         )
   return parser.parse_args(args)

# to retrieve command line arguments and assign to variables 
args = check_arg(sys.argv[1:])
infile = args.input
outfile = args.output

#establishing taxa list 
taxa_id_list = []

with open (infile, "r") as r:

    #getting only the taxa_id columns from READMe 
    taxa_array = np.loadtxt(r, usecols=1,skiprows=1,dtype=str)

#making the taxaID array into a list of taxaIDs  
for taxa_array in taxa_array: 
    taxa_id_list.append(str(taxa_array))

#to create the newick species tree, the tool ete4 was used 
#to install: pip install ete4 (in command line)
from ete4 import NCBITaxa


ncbi = NCBITaxa()
tree = ncbi.get_topology(taxa_id_list) #makes species tree with taxaIDs  

taxa_name_list = [int(node.name) for node in tree.leaves()] #making a list of taxaIDs in order  


trans_name = ncbi.get_taxid_translator(taxa_name_list) #translates IDs to species and stores IDs and species names in a dictionary 


#replace taxaID in tree node with name of species in the tree node 
for node in tree.traverse():
    taxa = int(node.name)
    if taxa in trans_name:
        node.name = trans_name[taxa].replace(" ", "_")

#outwritting the tree in newick formatting 
tree.write(parser=1, outfile="species_tree_2020.nw")

