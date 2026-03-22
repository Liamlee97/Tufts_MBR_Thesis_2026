"""
Figure 3.4 Pipeline - UMAP Visualization
----------------------------------------
This script isolates key marker genes for identified cellular populations
from the scVI integrated dataset, rendering custom normalized UMAP components
for publication.

Key Steps:
- Load the fully integrated `.h5ad` human dataset object.
- Select target marker genes (TP63, KRT5, NRCAM, ADH7).
- Generates high-resolution publication-themed UMAP visuals plotting those readouts.
"""

import scanpy as sc
import sys
import os

def main():
    # --------------------------------------------------------------------------
    # Configuration and Paths
    # --------------------------------------------------------------------------
    h5ad_file = "Concat_dataset_integrated.h5ad" 
    try:
        adata = sc.read_h5ad(h5ad_file)
    except FileNotFoundError:
        print(f"Error: {h5ad_file} missing from {os.path.abspath('.')}.")
        return

    # --------------------------------------------------------------------------
    # Plotting Formatting
    # --------------------------------------------------------------------------
    sc.settings.set_figure_params(figsize=(8, 6), dpi_save=600)

    genes = ['TP63', 'KRT5', 'NRCAM', 'ADH7']
    
    valid_genes = [g for g in genes if g in adata.var_names]

    if valid_genes:
        fig = sc.pl.umap(
            adata,
            color=valid_genes,
            ncols=2,
            s=3,
            legend_loc='right margin',
            legend_fontsize=16,
            layer='lognorm',
            show=False,
            return_fig=True
        )

        for ax in fig.axes:
            if ax.get_title() in valid_genes:
                ax.set_xlabel("UMAP1", fontsize=30)
                ax.set_ylabel("UMAP2", fontsize=30)
                ax.set_title(ax.get_title(), fontsize=30)

        out_name = "Figure 3.4_UMAP_MarkerGenes.tiff"
        fig.savefig(
            out_name,
            dpi=600,
            bbox_inches="tight"
        )
        print(f"Saved {out_name} successfully!")
    else:
        print("Required marker genes were absent from var_names list.")

if __name__ == "__main__":
    main()