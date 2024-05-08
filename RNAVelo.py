import numpy as np
import pandas as pd

from pyroe import load_fry
import scanpy as sc
import scvelo as scv
import anndata
from tqdm import tqdm
import pickle

import umap

import matplotlib.pyplot as plt
import os
import sys


srr_id = str(sys.argv[1])
n_cores = int(sys.argv[2])

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor="white")
sc.settings.n_jobs = 8


def fill_unmapped(name_df):
    '''Some ENTREZ IDs are not mapped to gene symbols. Fill them in or else downstream stuff complains'''
    na_i = name_df["GENE_SYMBOL"].isna()
    filler = [f"UNDEF{x}" for x in range(na_i.sum())]

    name_df.loc[name_df["GENE_SYMBOL"].isna(), "GENE_SYMBOL"] = filler
    
    return name_df


def number_duplicates(arr):
    '''Number duplicate entries so things don't complain'''
    cat_arr = pd.Series([""]*len(arr))
    dup_i = arr.duplicated()

    cat_arr[dup_i] = ".2"
    return arr.str.cat(cat_arr)


srr_sample = {
    "SRR26383970":"sample12",
    "SRR26383973":"sample15",
    "SRR26383988":"sample30",
    "SRR26383989":"sample31",
    "SRR26383990":"sample32",
    "SRR26383996":"sample38",
    "SRR26384000":"sample42",
    "SRR26383987":"sample29",
    "SRR26383992":"sample34",
    "SRR26383994":"sample36",
    "SRR26383999":"sample41"
}

md_df = pd.read_csv("./AlevinFry/aux_metadata_EL.csv")
smd_df = md_df[md_df["orig.ident"] == srr_sample[srr_id]]
smd_df["barcodes"] = smd_df["Unnamed: 0"].str.split("-").str[0]


# These genes come from the nature paper
ct_markers = {
    "Mesenchymal cells": ["DCN", "COL11A1", "PDGFRA","PDGFRB"],
    "Epithelial cells": ["EPCAM", "KRT8", "KRT10", "KRT18", "KRT19"],
    "Smooth muscle cells": ["ACTA2"],
    "Erythrocytes": ["HBB", "GYPA"],
    "Endothelial cells": ["CLDN5", "PECAM1", "CD34", "ESAM"],
    "Mast cells": ["TPSB2"],
    "Myeloid cells": ["LYZ", "CD14","C1QA", "CLEC10A"],
    "T/NK cells": ["CD2", "CD3D", "CD3E", "CD3G", 'CD8A', "CCL5"],
    "B/Plasma cells": ["JCHAIN", "CD79A"]
}

mgs = np.concatenate(list(ct_markers.values()))


# Load network data
mgs_ensg = list(pd.read_csv("./marker_genes.csv")["ENSEMBL"])
kegg_db = pd.read_csv("./kegg_pw_db.csv")

pw_ensmbl = list(kegg_db[["ENSEMBL", "Symbol"]].drop_duplicates()["ENSEMBL"])
pw_ensmbl.extend(mgs_ensg)
pw_ensmbl = np.unique(pw_ensmbl)


# load count matrix
adata = load_fry(f"/home/gridsan/elu/20.440Project/AlevinFry/quants/{srr_id}_quant_v3/af_quant/", output_format = "velocity")

# Load gene name mapping
name_df = pd.read_csv("./gene_name_mapping.csv")
name_df = fill_unmapped(name_df)

scv.pp.filter_and_normalize(adata, n_top_genes=5000)

# Streamline namings
merged = adata.var.merge(name_df, how="inner", left_index=True, right_on="ENSG_ID")
merged.index = merged["GENE_SYMBOL"]
merged.index = number_duplicates(pd.Series(merged.index))
merged.index.name = "gene_ids"
merged = merged.drop(["ENSG_ID", "GENE_SYMBOL"], axis=1)

adata.var = merged
adata.var.index = number_duplicates(pd.Series(adata.var.index))


# Filter out unimportant cells and the cells that authors filtered out
excluded_cts = ["T/NK cells", "B/Plasma cells", "Mast cells", "Erythrocytes", "Myeloid cells"]
smd_df = smd_df[~smd_df["active.cluster"].isin(excluded_cts)]

bc_mask = adata.obs["barcodes"].isin(smd_df["barcodes"])
adata = adata[bc_mask, :]


sc.pp.neighbors(adata)

# sc.pp.neighbors(adata)
scv.tl.umap(adata, random_state=1, min_dist=0.1)
umap_x = adata.obsm["X_umap"]


cts_df = adata.to_df()
avail_mgs = np.intersect1d(cts_df.columns, mgs)
print("Markers filtered out:",  ', '.join(np.setdiff1d(mgs, avail_mgs)))

# Record filtered-out genes
mg_df = cts_df[avail_mgs]
filtered_mgs = list(np.setdiff1d(mgs, avail_mgs))

# Remove filtered-out genes from ct_markers
for k, v in ct_markers.items():
    for i in range(len(filtered_mgs)):
        if filtered_mgs[i] in ct_markers[k]:
            ct_markers[k].remove(filtered_mgs[i])

# Louvain clustering
sc.tl.louvain(adata)

# Ranking and testing genes between groups
sc.tl.rank_genes_groups(adata, mask_var=avail_mgs, groupby='louvain')


# Full dynamics
full_dynamics = True
if full_dynamics:
    scv.tl.recover_dynamics(adata, n_jobs=int(n_cores//2), var_names="all")

# Expeditedt
scv.tl.velocity(adata, n_jobs=int(n_cores//2))

scv.tl.velocity_graph(adata)


fig, ax = plt.subplots(figsize=(10, 8))

u_cts = list(smd_df["active.cluster"].drop_duplicates())

for i in range(len(u_cts)):
    smd_cluster = smd_df[smd_df["active.cluster"] == u_cts[i]]["barcodes"]
    ct_pts = np.where(adata.obs["barcodes"].isin(smd_cluster))[0]
    ax.scatter(adata.obsm["X_umap"][ct_pts,0], adata.obsm["X_umap"][ct_pts,1], s=10, label=u_cts[i])
    
scv.tl.velocity_graph(adata, n_jobs=int(n_cores//2))
scv.pl.velocity_embedding_stream(adata, basis='umap', ax=ax, smooth=0.5, color="white", show=False)

lgnd = ax.legend(loc="lower center", ncols=3, bbox_to_anchor=(0.5, -0.25))
for handle in lgnd.legend_handles:
    handle.set_sizes([100])

plt.savefig(f"./figures/filt_{srr_id}_velo.jpg", dpi=300, bbox_inches="tight")
plt.show()


with open(f"./data/pickles/filt_{srr_id}.pkl", "wb") as f:
    pickle.dump(adata, f)