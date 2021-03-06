
# %%
import pandas as pd
import igraph
import os

path = "ds_2020_12_15_9_4_1/"

rin = []
ain = []

def r_a_INExtractor(path:str) -> tuple(list(), list()):
    '''
    Extracts atom and residue interaction networks from a collection of getcontacts RIN .tsv's in the specified path.
    '''
    # use the residuizer as a filter for the residue network
    residuizer = lambda x: str(x).split(":")[1]+"_"+str(x).split(":")[2]

    rin_res = []    # residue interaction network
    rin_atm = [] # atom interaction network
    fileNames = []
    for file in os.scandir(path):
        if (file.path.endswith("statcont_all.tsv") and file.is_file()):
            # append filename
            fileNames.append(file.name.split("_")[0])

            # atomic-level interaction
            df_a = pd.read_csv(file, sep='\t', skiprows=2,names=["frame", "interaction_type", "atom_1", "atom_2"])[["atom_1", "interaction_type", "atom_2"]]
            df_a.rename(columns={"atom_1":"atm1", "interaction_type":"iType", "atom_2":"atm2"}, inplace = True)

            rin_atm.append(df_a)

            # residue-level interaction
            df_r = pd.DataFrame([df_a["atm1"].apply(residuizer), df_a["iType"], df_a["atm2"].apply(residuizer)]).transpose()
            df_r.drop_duplicates()
            df_r.rename(columns={"atm1":"resA", "iType":"iType", "atm2":"resB"}, inplace = True)
            df_r.drop_duplicates()
            rin_res.append(df_r)
    
    return (rin_atm, rin_res)
    '''
    Extracts atom and residue interaction networks from a collection of getcontacts RIN .tsv's in the specified path.
    '''
    # use the residuizer as a filter for the residue network
    residuizer = lambda x: str(x).split(":")[1]+"_"+str(x).split(":")[2]

    rin_res = []    # residue interaction network
    rin_atm = [] # atom interaction network
    fileNames = []
    for file in os.scandir(path):
        if (file.path.endswith("statcont_all.tsv") and file.is_file()):
            # append filename
            fileNames.append(file.name.split("_")[0])

            def getAINdf(file: file):
                df = pd.read_csv(file, sep='\t', skiprows=2,names=["frame", "interaction_type", "atom_1", "atom_2"])[["atom_1", "interaction_type", "atom_2"]]
                df.rename(columns={"atom_1":"atm1", "interaction_type":"iType", "atom_2":"atm2"}, inplace = True)
                return df

            # atomic-level interaction
            df_a = pd.read_csv(file, sep='\t', skiprows=2,names=["frame", "interaction_type", "atom_1", "atom_2"])[["atom_1", "interaction_type", "atom_2"]]
            df_a.rename(columns={"atom_1":"atm1", "interaction_type":"iType", "atom_2":"atm2"}, inplace = True)

            def getAINam(df:pd.DataFrame) -> list():
                '''
                Returns adjacency matrix from an atomic-level interaction dataframe.
                '''
                a_el = df_a[['atm1','atm2']].values.tolist()  # edge list
                return igraph.Graph.TupleList(a_el).get_adjacency()

            rin_atm.append(df_a)

            # residue-level interaction
            df_r = pd.DataFrame([df_a["atm1"].apply(residuizer), df_a["iType"], df_a["atm2"].apply(residuizer)]).transpose()
            df_r.drop_duplicates()
            df_r.rename(columns={"atm1":"resA", "iType":"iType", "atm2":"resB"}, inplace = True)
            df_r.drop_duplicates()
            rin_res.append(df_r)
    
    return (rin_atm, rin_res)

(ain, rin) = r_a_IN_Extractor(path)

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

for df in rin:
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
#total_corr.style.background_gradient(cmap='coolwarm',axis=None)
total_corr.sort_index(ascending=False, axis=1).sort_index(ascending=False, axis=0)

total_corr.to_csv("total_correlation.csv")