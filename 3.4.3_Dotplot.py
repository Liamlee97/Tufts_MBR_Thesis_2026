"""
Figure 3.4 Pipeline - Differential Gene Expression and Visualization
------------------------------------------------------------------
This script visualizes integration analysis for Figure 3.4.
It utilizes the scVI-integrated latent representation to identify distinct 
functional Horizontal Basal Cell (HBC) and Respiratory Basal Cell 
(Resp BC) populations.

Key Steps:
- Load the scVI-integrated human dataset cluster map.
- Annotate specific high-resolution Leiden clusters as 'HBC' (cluster 18) 
  and 'Resp BC' (clusters 3, 8, 14, 26).
- Compute differential gene expression between HBC and Resp BC populations
  using a Wilcoxon rank-sum test on log-normalized read counts.
- Export the top differential expression results to a CSV log.
- Generate a publication-quality dotplot summarizing the top markers
  distinguishing HBC and Resp BC clusters.
"""

import scanpy as sc
import matplotlib.pyplot as plt
import os

def main():
    # --------------------------------------------------------------------------
    # Load the Integrated Dataset
    # --------------------------------------------------------------------------
    input_path = "Concat_dataset_integrated.h5ad"
    print(f"Loading integrated dataset from {input_path}...")
    
    try:
        adata = sc.read_h5ad(input_path)
    except FileNotFoundError:
        print(f"Error: Could not locate integrated data at {os.path.abspath(input_path)}")
        print("Ensure integration pipeline steps have been executed prior to running this script.")
        return

    # --------------------------------------------------------------------------
    # Annotate Cell Populations Based on Clusters
    # --------------------------------------------------------------------------
    print("Isolating HBC and Resp BC populations...")
    
    pop1_mask = (adata.obs['leiden_scVI_1.2'] == '18') 
    pop2_mask = adata.obs['leiden_scVI_1.2'].isin(['3', '8', '14', '26'])

    adata.obs['population'] = 'other'
    adata.obs.loc[pop1_mask, "population"] = "HBC"
    adata.obs.loc[pop2_mask, "population"] = "Resp BC"

    # --------------------------------------------------------------------------
    # Differential Gene Expression Analysis
    # --------------------------------------------------------------------------
    print("Computing differential gene expression (HBC vs Resp BC)...")
    
    sc.tl.rank_genes_groups(
        adata, 
        groupby='population', 
        groups=['HBC'], 
        reference='Resp BC', 
        method='wilcoxon',
        layer="lognorm"
    )

    results_df = sc.get.rank_genes_groups_df(adata, group='HBC')
    
    csv_path = "differential_leiden.csv"
    results_df.to_csv(csv_path, index=False)
    print(f"Differential gene expression results saved to {os.path.abspath(csv_path)}")

    # --------------------------------------------------------------------------
    # Generate Publication-Quality Visualizations
    # --------------------------------------------------------------------------
    print("Generating differential expression dotplot...")
    
    genes_top = results_df["names"].head(10).tolist()

    adata_plot = adata[adata.obs["population"].isin(["HBC", "Resp BC"])].copy()

    sc.settings.set_figure_params(figsize=(8, 6), dpi_save=600, fontsize=16)

    dp = sc.pl.dotplot(
        adata_plot,
        var_names=genes_top,
        groupby="population",
        categories_order=["HBC", "Resp BC"],
        layer='lognorm',
        dot_min=0.05,
        dot_max=0.8,
        figsize=(3, 5),          
        show=False,
        swap_axes=True,
        cmap="Reds",
        return_fig=True
    )

    fig = plt.gcf()
    fig.subplots_adjust(left=0.35, right=1.5)

    for ax in fig.axes:
        ax.tick_params(axis="y", labelsize=9)

    out_file = "Figure 3.4_HBC_vs_RespBC_dotplot.tiff"
    fig.savefig(out_file, dpi=600, bbox_inches="tight")
    print(f"Successfully generated and saved {out_file}.")

if __name__ == "__main__":
    main()