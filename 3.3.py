"""
Figure 3.3 Murine scRNA-seq analysis identifies candidate markers for HBC enrichment
------------------------------
This script generates the plots for Figure 3.3 using the public GSE245074 dataset.
It performs a sub-clustering specifically on the 'HBC' population to separate
true Horizontal Basal Cells (HBCs) from Respiratory Basal Cells (Resp BCs).

The script outputs:
1. UMAP highlighting the identified lineages (excluding immune cells).
2. UMAPs for key marker genes (Trp63, Krt5, Nrcam, Adh7).
3. A dotplot showing the top differentially expressed genes between HBC and Resp BC.

Requirements: GSE245074.h5ad
"""

import os
import scanpy as sc
import matplotlib.pyplot as plt

def main():
    # --------------------------------------------------------------------------
    # Load Data
    # --------------------------------------------------------------------------
    dataset_path = "GSE245074.h5ad"
    
    if not os.path.exists(dataset_path):
        print(f"Error: {dataset_path} not found.")
        print("Please ensure the publicly available GSE245074 processed data is in the current directory.")
        return

    print("Loading GSE245074...")
    adata = sc.read_h5ad(dataset_path)

    # --------------------------------------------------------------------------
    # Re-cluster the HBC population to separate HBCs from Resp BCs
    # --------------------------------------------------------------------------
    print("Sub-clustering HBC population...")
    
    if 'connectivities' not in adata.obsp or 'distances' not in adata.obsp:
        print("Neighbor graph not found. Computing neighbors...")
        sc.pp.neighbors(adata, n_neighbors=14, n_pcs=20)

    adata.obs['cell_labels'] = adata.obs['cell_labels'].astype('category')
    
    # Restrict Leiden clustering (resolution 0.2) only to the 'HBC' label
    try:
        sc.tl.leiden(adata, resolution=0.2, restrict_to=('cell_labels', ['HBC']), key_added='cell_labels_v2')
    except Exception as e:
        print(f"Error during Leiden clustering: {e}")
        print("Ensure that 'HBC' exists in adata.obs['cell_labels']")
        return

    # Map the resulting sub-clusters stringently based on established markers
    def map_hbc_subclusters(x):
        mapping = {
            "HBC,2": "HBC",
            "HBC,4": "Basal-Hillock-2",
            "HBC,0": "Resp BC",
            "HBC,1": "Resp BC",
            "HBC,3": "Resp BC",
            "HBC,5": "Resp BC",
            "HBC,6": "Resp BC",
            "HBC,7": "Resp BC"
        }
        return mapping.get(x, x)

    adata.obs['cell_labels_final'] = adata.obs['cell_labels_v2'].map(map_hbc_subclusters).astype("category")

    # --------------------------------------------------------------------------
    # Filter Immune Cells for Visualization
    # --------------------------------------------------------------------------
    print("Filtering immune cells...")
    immune_labels = [
        "B Cell", "pre-B Cell", "Plasma Cell", "T cell (CD4)", "T cell (CD8)",
        "Naive T cell (CD4)", "Naive T cell (CD8)", "gd T cell", "TH2", "TH17",
        "Treg", "NKT", "NK", "NK-CD27", "ILC1", "ILC2", "ILC3", "Monocyte",
        "Macrophage", "DC", "pDC", "Neutrophil", "pre-Neutrophil", "Eosinophil",
        "Basophil", "Mast Cell", "CMP", "IEL", "Osteoclast"
    ]

    adata_no_immune = adata[~adata.obs["cell_labels_final"].isin(immune_labels)].copy()
    
    # --------------------------------------------------------------------------
    # Plot Generation
    # --------------------------------------------------------------------------
    sc.settings.set_figure_params(figsize=(8, 6), dpi_save=600, fontsize=16)
    
    print("Generating Figure 3.3a: Cell Type UMAP...")
    fig_umap = sc.pl.umap(
        adata_no_immune, 
        color=['cell_labels_final'],
        ncols=1,
        legend_loc='right margin',
        legend_fontsize=16,
        show=False,
        return_fig=True
    )
    
    for ax in fig_umap.axes:
        ax.set_xlabel("UMAP1", fontsize=16)
        ax.set_ylabel("UMAP2", fontsize=16)
        ax.set_title("GSE245074 (Immune Subset Filtered)", fontsize=20)

    fig_umap.savefig("Figure 3.3_UMAP_CellTypes.tiff", dpi=600, bbox_inches="tight")
    print("Saved Figure 3.3_UMAP_CellTypes.tiff")

    print("Generating Figure 3.3b: Marker Gene UMAPs...")
    genes = ['Trp63', 'Krt5', 'Nrcam', 'Adh7']
    
    valid_genes = [g for g in genes if g in adata_no_immune.var_names]
    
    if valid_genes:
        fig_genes = sc.pl.umap(
            adata_no_immune, 
            color=valid_genes,
            ncols=2,
            legend_loc='right margin',
            legend_fontsize=16,
            show=False,
            return_fig=True
        )

        for ax in fig_genes.axes:
            if ax.get_title() in valid_genes:
                ax.set_xlabel("UMAP1", fontsize=30)
                ax.set_ylabel("UMAP2", fontsize=30)
                ax.set_title(ax.get_title(), fontsize=30)

        fig_genes.savefig("Figure 3.3_UMAP_MarkerGenes.tiff", dpi=600, bbox_inches="tight")
        print("Saved Figure 3.3_UMAP_MarkerGenes.tiff")
    else:
        print("Warning: None of the target marker genes were found in the dataset.")

    print("Generating Figure 3.3c: DEG Dotplot...")
    
    adata_sub = adata_no_immune[adata_no_immune.obs["cell_labels_final"].isin(["HBC", "Resp BC"])].copy()
    
    if "HBC" in adata_sub.obs["cell_labels_final"].values and "Resp BC" in adata_sub.obs["cell_labels_final"].values:
        print("Calculating differentially expressed genes (HBC vs Resp BC)...")
        sc.tl.rank_genes_groups(
            adata_sub,
            groupby="cell_labels_final",
            groups=["HBC"],
            reference="Resp BC",
            method="wilcoxon"
        )

        deg = sc.get.rank_genes_groups_df(adata_sub, group="HBC")
        genes_top = deg["names"].head(10).tolist()
        
        print(f"Top 10 markers for HBC: {genes_top}")

        group_order = ["HBC", "Resp BC"]

        dp = sc.pl.dotplot(
            adata_sub,
            var_names=genes_top,
            groupby="cell_labels_final",
            categories_order=group_order,
            dot_min=0.05,
            dot_max=0.8,
            figsize=(3, 5),
            show=False,
            swap_axes=True,
            cmap="Reds",
            return_fig=True
        )

        dotplot_fig = plt.gcf()
        dotplot_fig.subplots_adjust(left=0.35, right=1.5)

        for ax in dotplot_fig.axes:
            ax.tick_params(axis="y", labelsize=9)

        dotplot_fig.savefig("Figure 3.3_Dotplot_HBC_vs_RespBC.tiff", dpi=600, bbox_inches="tight")
        print("Saved Figure 3.3_Dotplot_HBC_vs_RespBC.tiff")
    else:
        print("Warning: Could not find both 'HBC' and 'Resp BC' for DEG computation.")
    
    print("Figure 3.3 generation complete.")

if __name__ == '__main__':
    main()
