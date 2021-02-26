# ring_analyzer.py

# %% 
# Imports
import pandas as pd
import numpy as np
import igraph
pd.set_option("display.precision", 3)
# %%
# get protein structures
# Libraries
import pypdb       # For searching and downloading PDB files
import datetime     # for organizing file output
import Bio      # Bio for downloading PDBs, printing protein lengths in residues
from Bio.PDB import PDBList
from Bio.PDB import PDBParser

# %%
query = input("Supply a query (term, accession number, etc.) :")

query = "nuclear receptor"
# Tag the time right when the query is entered

now = datetime.datetime.now()
def now_dir_ts():
    now_ts = str(now.year)+"_"+str(now.month)+"_"+str(now.day)+"_"+str(now.hour)+"_"+str(now.minute)+"_"+str(now.second)
    return now_ts

now = now_dir_ts()
PDB_dl_dir = "ds_"+now

search_dict = pypdb.Query(query)     # create a dictionary containing search information
found = search_dict.search(search_dict)[:500]      # create a list of these PDBs by searching RCSB

# create a list with the information and the metadata
metadata = []

for proteins in found:  # for items in # for the items in the list,
    metadata.append(pypdb.describe_pdb(proteins))  # append the dictionary 

# Save the metadata list as a CSV file
dfm = pd.DataFrame(metadata) # convert to a Pandas DF
dfm.to_csv('metadata_'+now+'.csv')      # save as a CSV file

# %%
parser = PDBParser()       # create a parser
pdbs = list()
pdbl = PDBList()

# Download all PDB structures in the previous list if they aren't there
for id in found:
    pdbl.retrieve_pdb_file(pdb_code=id, file_format='pdb', pdir=PDB_dl_dir)   # Retrieve in PDB format, put in directory 'PDB'

# Finished, print "Downloading ... finished!"
print('\n#############~DOWNLOADING COMPLETE~#############\n')
# %%
# convert pdb*.ent to *.pdb
for file in os.scandir(PDB_dl_dir):
    if (file.path.endswith(".ent") and file.is_file()):
        newfn = file.name.replace("pdb","").replace(".ent",".pdb")
        os.rename(file, PDB_dl_dir+"/"+newfn)
