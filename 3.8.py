"""
Figure 3.8. LIANA+ analysis identifies putative signaling pathways between HBCs and ionocytes.
--------------------------------------------------
This script performs a rank-aggregated cell-cell communication analysis
using LIANA to identify ligand-receptor interactions between Ionocyte 
and HBC populations from single-cell RNA-seq data.

Outputs: IonocytetoHBC.svg, HBCtoIonocyte.svg
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc
import liana as li
import decoupler as dc
import anndata as ad
from liana.method import cellphonedb

# -----------------------------------------------------------------------------
# Figure Settings & Data Loading
# -----------------------------------------------------------------------------
sc.settings.set_figure_params(dpi=300, frameon=False)
sc.set_figure_params(dpi=300, facecolor="white")
sc.set_figure_params(figsize=(5, 5))

file_path = 'GSE245074.h5ad'
print(f"Loading {file_path}...")
try:
    adata = sc.read_h5ad(file_path)
except FileNotFoundError:
    print(f"ERROR: {file_path} not found.")
    print("Please download the public dataset from GEO: GSE245074")
    import sys; sys.exit(1)

print("\nCell Label Frequencies:")
print(adata.obs["cell_labels"].value_counts().head(20))

# -----------------------------------------------------------------------------
# LIANA Setup & Ortholog Mapping (Human -> Mouse)
# -----------------------------------------------------------------------------
resource = li.rs.select_resource('consensus')

print("\nFetching HCOP Orthologs...")
map_df = li.rs.get_hcop_orthologs(url='https://ftp.ebi.ac.uk/pub/databases/genenames/hcop/human_mouse_hcop_fifteen_column.txt.gz',
                                  columns=['human_symbol', 'mouse_symbol'],
                                  min_evidence=3)

map_df = map_df.rename(columns={'human_symbol':'source', 'mouse_symbol':'target'})

mouse = li.rs.translate_resource(resource,
                                 map_df=map_df,
                                 columns=['ligand', 'receptor'],
                                 replace=True,
                                 one_to_many=1)

# -----------------------------------------------------------------------------
# Rank Aggregation & Plotting
# -----------------------------------------------------------------------------
print("\nRunning LIANA Rank Aggregate...")
li.mt.rank_aggregate(
    adata,
    groupby="cell_labels",
    resource_name='consensus',
    expr_prop=0.1,
    use_raw=False,
    verbose=True, 
    resource=mouse
)

plot_ionocyte_hbc = li.pl.dotplot(adata,
              colour='magnitude_rank',
              size='specificity_rank',
              inverse_size=True,
              inverse_colour=True,
              source_labels=['Ionocyte'],
              target_labels=['HBC'],
              top_n=30,
              orderby='magnitude_rank',
              orderby_ascending=True,
              figure_size=(5, 4),
              filter_fun=lambda x: x['specificity_rank'] <= 0.02,
             )
plot_ionocyte_hbc.save(filename="IonocytetoHBC.svg")

plot_hbc_ionocyte = li.pl.dotplot(adata,
              colour='magnitude_rank',
              size='specificity_rank',
              inverse_size=True,
              inverse_colour=True,
              source_labels=['HBC'],
              target_labels=['Ionocyte'],
              top_n=30,
              orderby='magnitude_rank',
              orderby_ascending=True,
              figure_size=(5, 4),
              filter_fun=lambda x: x['specificity_rank'] <= 0.02,
             )
plot_hbc_ionocyte.save(filename="HBCtoIonocyte.svg")

print("\nAnalysis complete. Figures saved as IonocytetoHBC.svg and HBCtoIonocyte.svg")
