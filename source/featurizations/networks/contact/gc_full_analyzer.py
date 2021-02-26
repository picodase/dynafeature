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

# %%
# run getcontacts on each pdb file in the folder
!qsub qsub_get_structnets.sh

#%%
import os

path = PDB_dl_dir

residuizer = lambda x: str(x).split(":")[1]+"_"+str(x).split(":")[2]

rin_res = []
rin_atomic = []
fileNames = []
for file in os.scandir(path):
    if (file.path.endswith("statcont_all.tsv") and file.is_file()):
        # append filename
        fileNames.append(file.name.split("_")[0])

        # atomic-level interaction
        df_a = pd.read_csv(file, sep='\t', skiprows=2,names=["frame", "interaction_type", "atom_1", "atom_2"])[["atom_1", "interaction_type", "atom_2"]]
        df_a.rename(columns={"atom_1":"atm1", "interaction_type":"iType", "atom_2":"atm2"}, inplace = True)

        rin_atomic.append(df_a)

        # residue-level interaction
        df_r = pd.DataFrame([df_a["atm1"].apply(residuizer), df_a["iType"], df_a["atm2"].apply(residuizer)]).transpose()
        df_r.drop_duplicates()
        df_r.rename(columns={"atm1":"resA", "iType":"iType", "atm2":"resB"}, inplace = True)
        df_r.drop_duplicates()
        rin_res.append(df_r)

# %%
# partitioning scheme for each network type

# types are hp,sb,pc,ps,ts,vdw,hb

tot = {"hp","sb","pc","ps","ts","vdw","hb"}
hbond = {"hb"}
vdw = {"vdw"}
sb = {"sb"}
πs = {"ps"}
πc = {"pc"}
ts = {"ts"}
hp = {"hp"}

dicts = [
    tot,
    hbond,
    vdw,
    sb,
    πs,
    πc,
    ts,
    hp
]

# %%
# create lists to hold each set of values for the specified interaction types
interac_cts = []
edgelists = []
graphs = []
#d_norm = []
thr_motifs = []
four_motifs = []

def constrInteracCts(d:dict):
    return df[df.iType.isin(d)][['resA', 'resB']].values.tolist()

for df in rin_res:
    # interaction counts
    interac_cts.append(df.iType.value_counts())

    # interaction types
    i = {}
    i["tot"] = df[['resA','resB']].values.tolist()
    i["hbd"] = df[df.iType.isin(hbond)][['resA','resB']].values.tolist()
    i["vdw"] = df[df.iType.isin(vdw)][['resA','resB']].values.tolist()
    i["ππ"] = df[df.iType.isin(πs)][['resA','resB']].values.tolist()
    i["sb"] = df[df.iType.isin(sb)][['resA','resB']].values.tolist()
    i["ps"] = df[df.iType.isin(πs)][['resA','resB']].values.tolist()
    i["pc"] = df[df.iType.isin(πc)][['resA','resB']].values.tolist()
    i["ts"] = df[df.iType.isin(ts)][['resA','resB']].values.tolist()
    i["hp"] = df[df.iType.isin(hp)][['resA','resB']].values.tolist()
    edgelists.append(i)

    # graphs
    g = {}
    g["tot"] = igraph.Graph.TupleList(i["tot"])
    g["hbd"] = igraph.Graph.TupleList(i["hbd"], directed=True)
    g["vdw"] = igraph.Graph.TupleList(i["vdw"])
    g["ππ"] = igraph.Graph.TupleList(i["ππ"])
    g["sb"] = igraph.Graph.TupleList(i["sb"])
    g["ps"] = igraph.Graph.TupleList(i["ps"])
    g["pc"] = igraph.Graph.TupleList(i["pc"])
    g["ts"] = igraph.Graph.TupleList(i["ts"])
    g["hp"] = igraph.Graph.TupleList(i["hp"])
    graphs.append(g)

    # motifs
    # calculate motifs for each graph

    # three-motifs
    t_m = {}
    t_m["tot"] = g["tot"].motifs_randesu()
    t_m["hbd"] = g["hbd"].motifs_randesu()
    t_m["vdw"] = g["vdw"].motifs_randesu()
    t_m["ππ"] = g["ππ"].motifs_randesu()
    t_m["sb"] = g["sb"].motifs_randesu()    
    t_m["ps"] = g["ps"].motifs_randesu()    
    t_m["pc"] = g["pc"].motifs_randesu()    
    t_m["ts"] = g["ts"].motifs_randesu()    
    t_m["hp"] = g["hp"].motifs_randesu()    
    thr_motifs.append(t_m)

    # four-motifs
    f_m = {}
    f_m["tot"] = g["tot"].motifs_randesu(size=4)
    f_m["hbd"] = g["hbd"].motifs_randesu(size=4)
    f_m["vdw"] = g["vdw"].motifs_randesu(size=4)
    f_m["ππ"] = g["ππ"].motifs_randesu(size=4)
    f_m["sb"] = g["sb"].motifs_randesu(size=4)    
    f_m["ps"] = g["ps"].motifs_randesu(size=4)    
    f_m["pc"] = g["pc"].motifs_randesu(size=4)    
    f_m["ts"] = g["ts"].motifs_randesu(size=4)    
    f_m["hp"] = g["hp"].motifs_randesu(size=4)  
    four_motifs.append(f_m)

# %%
import numpy as np
import sklearn
from sklearn.metrics.pairwise import cosine_similarity
import jinja2

def dfCosSim(n1: np.ndarray, n2: np.ndarray):
  return cosine_similarity(n1.reshape(1,-1), n2.reshape(1,-1))
# %%
# residue interaction counts
rin_interact_data = pd.DataFrame(interac_cts, index=fileNames).T
rn_interact_cts_corr = rin_interact_data.corr(dfCosSim)
rn_interact_cts_corr.style.background_gradient(cmap='coolwarm')
# %%
#Cosine Correlation Matrices for 3 motif

#thr_motif_tot = [d['tot'] for d in thr_motifs]
thr_motif_hbd = [d['hbd'] for d in thr_motifs]
thr_motif_vdw = [d['vdw'] for d in thr_motifs]
thr_motif_pipi = [d['ππ'] for d in thr_motifs]
thr_motif_sb = [d['sb'] for d in thr_motifs]
thr_motif_ps = [d['ps'] for d in thr_motifs]
thr_motif_pc = [d['pc'] for d in thr_motifs]
thr_motif_ts = [d['ts'] for d in thr_motifs]
thr_motif_hp = [d['hp'] for d in thr_motifs]

#rn_thr_motif_tot_data = pd.DataFrame(thr_motif_tot, index=fileNames).T
rn_thr_motif_hbd_data = pd.DataFrame(thr_motif_hbd, index=fileNames).T
rn_thr_motif_vdw_data = pd.DataFrame(thr_motif_vdw, index=fileNames).T
rn_thr_motif_pipi_data = pd.DataFrame(thr_motif_pipi, index=fileNames).T
rn_thr_motif_sb_data = pd.DataFrame(thr_motif_sb, index=fileNames).T
rn_thr_motif_ps_data = pd.DataFrame(thr_motif_ps, index=fileNames).T
rn_thr_motif_pc_data = pd.DataFrame(thr_motif_pc, index=fileNames).T
rn_thr_motif_ts_data = pd.DataFrame(thr_motif_ts, index=fileNames).T
rn_thr_motif_hp_data = pd.DataFrame(thr_motif_hp, index=fileNames).T

#rn_thr_motif_tot_corr = rn_thr_motif_tot_data.corr(dfCosSim)
rn_thr_motif_hbd_corr = rn_thr_motif_hbd_data.corr(dfCosSim)
rn_thr_motif_vdw_corr = rn_thr_motif_vdw_data.corr(dfCosSim)
rn_thr_motif_pipi_corr = rn_thr_motif_pipi_data.corr(dfCosSim)
rn_thr_motif_sb_corr = rn_thr_motif_sb_data.corr(dfCosSim)
rn_thr_motif_ps_corr = rn_thr_motif_ps_data.corr(dfCosSim)
rn_thr_motif_pc_corr = rn_thr_motif_pc_data.corr(dfCosSim)
rn_thr_motif_ts_corr = rn_thr_motif_ts_data.corr(dfCosSim)
rn_thr_motif_hp_corr = rn_thr_motif_hp_data.corr(dfCosSim)

# %%
#Cosine Correlation Matrices for 3 motif

#four_motif_tot = [d['tot'] for d in four_motifs]
four_motif_hbd = [d['hbd'] for d in four_motifs]
four_motif_vdw = [d['vdw'] for d in four_motifs]
four_motif_pipi = [d['ππ'] for d in four_motifs]
four_motif_sb = [d['sb'] for d in four_motifs]
four_motif_ps = [d['ps'] for d in four_motifs]
four_motif_pc = [d['pc'] for d in four_motifs]
four_motif_ts = [d['ts'] for d in four_motifs]
four_motif_hp = [d['hp'] for d in four_motifs]

#rn_four_motif_tot_data = pd.DataFrame(four_motif_tot, index=fileNames).T
rn_four_motif_hbd_data = pd.DataFrame(four_motif_hbd, index=fileNames).T
rn_four_motif_vdw_data = pd.DataFrame(four_motif_vdw, index=fileNames).T
rn_four_motif_pipi_data = pd.DataFrame(four_motif_pipi, index=fileNames).T
rn_four_motif_sb_data = pd.DataFrame(four_motif_sb, index=fileNames).T
rn_four_motif_ps_data = pd.DataFrame(four_motif_ps, index=fileNames).T
rn_four_motif_pc_data = pd.DataFrame(four_motif_pc, index=fileNames).T
rn_four_motif_ts_data = pd.DataFrame(four_motif_ts, index=fileNames).T
rn_four_motif_hp_data = pd.DataFrame(four_motif_hp, index=fileNames).T

#rn_four_motif_tot_corr = rn_four_motif_tot_data.corr(dfCosSim)
rn_four_motif_hbd_corr = rn_four_motif_hbd_data.corr(dfCosSim)
rn_four_motif_vdw_corr = rn_four_motif_vdw_data.corr(dfCosSim)
rn_four_motif_pipi_corr = rn_four_motif_pipi_data.corr(dfCosSim)
rn_four_motif_sb_corr = rn_four_motif_sb_data.corr(dfCosSim)
rn_four_motif_ps_corr = rn_four_motif_ps_data.corr(dfCosSim)
rn_four_motif_pc_corr = rn_four_motif_pc_data.corr(dfCosSim)
rn_four_motif_ts_corr = rn_four_motif_ts_data.corr(dfCosSim)
rn_four_motif_hp_corr = rn_four_motif_hp_data.corr(dfCosSim)

# %%
total_corr = (rn_interact_cts_corr + rn_thr_motif_hbd_corr + rn_thr_motif_vdw_corr + rn_thr_motif_ps_corr + rn_thr_motif_sb_corr + rn_thr_motif_ps_corr + rn_thr_motif_pc_corr + rn_thr_motif_ts_corr + rn_thr_motif_hp_corr + rn_four_motif_hbd_corr + rn_four_motif_vdw_corr + rn_four_motif_pipi_corr + rn_four_motif_sb_corr + rn_four_motif_ps_corr + rn_four_motif_pc_corr + rn_four_motif_ts_corr + rn_four_motif_hp_corr) / 17

# %%
renamed_total_corr.style.background_gradient(cmap='coolwarm',axis=None)

renamed_total_corr.to_csv("total_correlation.csv")