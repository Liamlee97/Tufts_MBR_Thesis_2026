"""
Figure 3.4 Pipeline - Integration using scVI
--------------------------------------------
This script utilizes the deep generative modeling framework scVI 
(single-cell Variational Inference) to computationally remove batch effects 
across the combined human datasets for Figure 3.4.

Key Steps:
- Select Highly Variable Genes (HVG) using scVI's Poisson flavor.
- Setup the scVI model parameters (accounting for batch via 'dataset' subset).
- Train the scVI generative model on raw counts for 500 epochs.
- Extract the latent biological representation representing the integrated data.
- Compute new neighborhood graphs and UMAP embedding on the integrated space.
- Perform Leiden clustering (resolution 1.2) for downstream differential analysis.
"""

import os
import scanpy as sc
import scvi
import pandas as pd
import numpy as np
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def main():
    # --------------------------------------------------------------------------
    # Configuration and Paths
    # --------------------------------------------------------------------------
    input_path = "Concat_dataset.h5ad"
    output_dir = "."
    
    logging.info(f"Reading concatenated data from {input_path}...")
    try:
        adata = sc.read_h5ad(input_path)
    except FileNotFoundError:
        logging.error(f"Input file {input_path} not found in {os.path.abspath(output_dir)}.")
        return
    
    # --------------------------------------------------------------------------
    # Feature Selection (Highly Variable Genes)
    # --------------------------------------------------------------------------
    batch_key = "dataset"
    logging.info(f"Available features/metadata columns: {adata.obs.columns.tolist()}")
    
    logging.info("Selecting Highly Variable Genes (HVG) via Seurat v3 flavor...")
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=3000,
        subset=False,
        layer="counts",
        flavor="seurat_v3",
        batch_key=batch_key
    )
    
    logging.info("Selecting HVG via scVI Poisson flavor...")
    df_poisson = scvi.data.poisson_gene_selection(
        adata, n_top_genes=3000, batch_key=batch_key, inplace=False
    )
    
    is_hvg = df_poisson.highly_variable
    adata.varm['df_poisson'] = df_poisson
    
    adata_query = adata[:, is_hvg].copy()
    logging.info(f"Subsetting to {adata_query.shape[1]} top Highly Variable Genes.")
    
    # --------------------------------------------------------------------------
    # Setup and Train the scVI Integration Model
    # --------------------------------------------------------------------------
    logging.info("Setting up scVI AnnData Registration...")
    scvi.model.SCVI.setup_anndata(
        adata_query,
        layer="counts",
        categorical_covariate_keys=[batch_key],
        continuous_covariate_keys=["pct_counts_mito"]
    )
    
    logging.info("Training the scVI generative model on negative binomial likelihood...")
    model = scvi.model.SCVI(adata_query, gene_likelihood="nb")
    
    train_kwargs = dict(
        early_stopping=True,
        early_stopping_patience=20,
        enable_model_summary=True,
        max_epochs=500
    )
    
    model.train(**train_kwargs)
    
    model_save_path = os.path.join(output_dir, "scvi_model")
    logging.info(f"Saving trained model to {os.path.abspath(model_save_path)}...")
    model.save(model_save_path, overwrite=True)
    
    # --------------------------------------------------------------------------
    # Extract Latent Representation and Cluster
    # --------------------------------------------------------------------------
    logging.info("Extracting latent representation array from trained model...")
    latent = model.get_latent_representation()
    adata.obsm["X_scVI"] = latent
    
    logging.info("Computing integration neighborhood graph and UMAP...")
    sc.pp.neighbors(adata, use_rep="X_scVI")
    sc.tl.umap(adata, min_dist=0.5)
    
    logging.info("Running high-resolution Leiden clustering (res=1.2)...")
    sc.tl.leiden(adata, key_added="leiden_scVI_1.2", resolution=1.2)
    
    # --------------------------------------------------------------------------
    # Save Final Integrated Object
    # --------------------------------------------------------------------------
    final_output_path = os.path.join(output_dir, 'Concat_dataset_integrated.h5ad')
    logging.info(f"Saving integrated dataset locally to {os.path.abspath(final_output_path)}...")
    adata.write(final_output_path)
    
    logging.info("Integration Pipeline Complete.")

if __name__ == "__main__":
    main()
