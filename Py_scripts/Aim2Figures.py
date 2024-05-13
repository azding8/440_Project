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

import seaborn as sns


sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor="white")


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


srr_id = "SRR26383973"

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

srr_df = pd.DataFrame(pd.Series(srr_sample)).reset_index().rename(columns={"index": "SRA_ID", 0: "orig.ident"})

md_df = md_df.merge(srr_df, on="orig.ident")
sample_df = md_df[["SRA_ID", "Anatomic annotation", "orig.ident", "Major.Class"]].drop_duplicates()


# load count matrix
adata = load_fry(f"/home/gridsan/elu/20.440Project/AlevinFry/quants/{srr_id}_quant_v3/af_quant/", output_format = "velocity")

# Load gene name mapping
name_df = pd.read_csv("./gene_name_mapping.csv")
name_df = fill_unmapped(name_df)

scv.pp.filter_and_normalize(adata, n_top_genes=2000)

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


u_cts = list(smd_df["active.cluster"].drop_duplicates())


fig, ax = plt.subplots(figsize=(10, 8))

u_cts = list(smd_df["active.cluster"].drop_duplicates())

for i in range(len(u_cts)):
    smd_cluster = smd_df[smd_df["active.cluster"] == u_cts[i]]["barcodes"]
    ct_pts = np.where(adata.obs["barcodes"].isin(smd_cluster))[0]
    ax.scatter(adata.obsm["X_umap"][ct_pts,0], adata.obsm["X_umap"][ct_pts,1], s=10, label=u_cts[i])

lgnd = ax.legend(loc="lower center", ncols=3, bbox_to_anchor=(0.5, -0.25))
for handle in lgnd.legend_handles:
    handle.set_sizes([100])

plt.show()


# Full dynamics
full_dynamics = True
if full_dynamics:
    scv.tl.recover_dynamics(adata, n_jobs=4, var_names="all")

# Expeditedt
scv.tl.velocity(adata, n_jobs=4)

scv.tl.velocity_graph(adata)


fig, ax = plt.subplots(figsize=(10, 8))

u_cts = list(smd_df["active.cluster"].drop_duplicates())

for i in range(len(u_cts)):
    smd_cluster = smd_df[smd_df["active.cluster"] == u_cts[i]]["barcodes"]
    ct_pts = np.where(adata.obs["barcodes"].isin(smd_cluster))[0]
    ax.scatter(adata.obsm["X_umap"][ct_pts,0], adata.obsm["X_umap"][ct_pts,1], s=3, label=u_cts[i])

scv.tl.velocity_graph(adata, n_jobs=4)
scv.pl.velocity_embedding_stream(adata, basis='umap', ax=ax, smooth=0.5, color="white", show=False)

lgnd = ax.legend(loc="lower center", ncols=3, bbox_to_anchor=(0.5, -0.25))
for handle in lgnd.legend_handles:
    handle.set_sizes([100])

plt.show()


pkls = os.listdir("./data/pickles/")
sra_ids = [x.split("_")[1].split(".")[0] for x in pkls]

var_dfs = dict()
for i in tqdm(range(len(pkls))):
    adata = pd.read_pickle(f"./data/pickles/{pkls[i]}")
    var_dfs[sra_ids[i]] = adata.var

df_list = list(var_dfs.values())
for i in range(len(sra_ids)):
    df = df_list[i]
    df["SRA_ID"] = sra_ids[i]
var_df = pd.concat(df_list)


hq_df = var_df[(var_df["fit_likelihood"] > 0.1) & var_df["velocity_genes"]].dropna().reset_index()

hq_df = hq_df.merge(sample_df)


splice_dict = {
    "gene_ids": [],
    "Unaffected ovary": [],
    "Endometrioma": [],
    "Eutopic Endometrium": []
}

gb_obj = hq_df.groupby("gene_ids")
for g, df in gb_obj:
    if len(df["Major.Class"].drop_duplicates()) > 1:
        splice_dict["gene_ids"].append(g)
        u_mcs = list(df["Major.Class"].drop_duplicates())
        for mc in ["Unaffected ovary", "Endometrioma", "Eutopic Endometrium"]:
            sub_df = df[df["Major.Class"] == mc]
            splice_dict[mc].append(sub_df["fit_beta"].mean())


splice_df = pd.DataFrame(splice_dict).dropna(subset="Unaffected ovary")
splice_df.index = splice_df["gene_ids"]
splice_df = splice_df.drop(columns="gene_ids")

splice_df = splice_df.replace({np.nan: 0})

gb_size = hq_df.groupby("gene_ids").size()

velo_genes = var_df[var_df["velocity_genes"]].dropna()
# velo_genes.groupby("gene_ids").size().value_counts()

velo_genes = velo_genes.reset_index()

velo_genes = velo_genes[velo_genes["fit_likelihood"] > 0.1]

# velo_genes = velo_genes[velo_genes["fit_pval_steady"] < 0.05]

velo_genes = velo_genes.merge(sample_df)
velo_genes["CC"] = velo_genes["Major.Class"] == "Unaffected ovary"


gb_obj = velo_genes.groupby("gene_ids")

n_velo_dict = {
    "Unaffected ovary": [],
    "Endometrioma": [],
    "Eutopic Endometrium": []
}

for g, df in gb_obj:
    n_velo_dict["Unaffected ovary"].append(len(df[df["Major.Class"] == "Unaffected ovary"]))
    n_velo_dict["Endometrioma"].append(len(df[df["Major.Class"] == "Endometrioma"]))
    n_velo_dict["Eutopic Endometrium"].append(len(df[df["Major.Class"] == "Eutopic Endometrium"]))


n_velo_df = pd.DataFrame(n_velo_dict)
n_velo_df["Unaffected ovary"] = n_velo_df["Unaffected ovary"]/3
n_velo_df["Endometrioma"] = n_velo_df["Endometrioma"]/4
n_velo_df["Eutopic Endometrium"] = n_velo_df["Eutopic Endometrium"]/4
n_velo_df.index = list(gb_obj.groups.keys())


gb_obj = velo_genes.groupby("gene_ids")

splice_dict = {
    "Unaffected ovary": [],
    "Endometrioma": [],
    "Eutopic Endometrium": []
}

for g, df in gb_obj:
    splice_dict["Unaffected ovary"].append(df[df["Major.Class"] == "Unaffected ovary"]["fit_beta"].mean())
    splice_dict["Endometrioma"].append(df[df["Major.Class"] == "Endometrioma"]["fit_beta"].mean())
    splice_dict["Eutopic Endometrium"].append(df[df["Major.Class"] == "Eutopic Endometrium"]["fit_beta"].mean())


splice_df = pd.DataFrame(splice_dict)
splice_df.index = list(gb_obj.groups.keys())


gb_obj = velo_genes.groupby("gene_ids")

size_dict = {
    "Unaffected ovary": [],
    "Endometrioma": [],
    "Eutopic Endometrium": []
}

for g, df in gb_obj:
    size_dict["Unaffected ovary"].append(df[df["Major.Class"] == "Unaffected ovary"].shape[0])
    size_dict["Endometrioma"].append(df[df["Major.Class"] == "Endometrioma"].shape[0])
    size_dict["Eutopic Endometrium"].append(df[df["Major.Class"] == "Eutopic Endometrium"].shape[0])

size_df = pd.DataFrame(size_dict)

n_velo_df = n_velo_df[n_velo_df.sum(axis=1) > 0.25]

sns.clustermap(splice_df, cmap="Reds", standard_scale=0, figsize=(5, 5))

plt.savefig("./figures/diff_splice.jpg", dpi=300, bbox_inches="tight")
plt.show()
