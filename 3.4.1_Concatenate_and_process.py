"""
Figure 3.4 Pipeline - Dataset Concatenation and Preprocessing
-------------------------------------------------------------
This script represents the initial integration pipeline for Figure 3.4.
It loads distinct human single-cell olfactory datasets, performs initial 
quality control filtering, calculates essential metrics (e.g., mitochondrial 
percentage), and merges them into a single concatenated raw object for 
downstream integration.

Key Steps:
- Load raw .h5ad files for separate datasets.
- Filter low-quality cells (< 200 genes) to optimize memory.
- Concatenate datasets, prefixing barcodes to preserve uniqueness.
- Calculate standard QC metrics (ribosomal/mitochondrial content).
- Filter cells with high mitochondrial content (> 20%).
- Generate log-normalized counts and save the processed checkpoint.
"""

import os
import scanpy as sc
import numpy as np
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def main():
    # --------------------------------------------------------------------------
    # Configuration and Paths
    # --------------------------------------------------------------------------
    human_raw_files = {
        "GSE139522": "GSE139522_raw.h5ad",
        "GSE184117": "GSE184117_normosmia_presbyosmia_pooled_raw.h5ad",
        "GSE201620": "GSE201620_raw.h5ad"
    }

    output_dir = "."

    # --------------------------------------------------------------------------
    # Load and Concatenate
    # --------------------------------------------------------------------------
    logging.info("Loading and concatenating 3 human raw datasets...")
    adatas = []
    
    for name, path in human_raw_files.items():
        logging.info(f"Loading {name} from {path}...")
        try:
            adata = sc.read_h5ad(os.path.join(output_dir, path))
        except FileNotFoundError:
            logging.error(f"Required raw file {path} not found in {os.path.abspath(output_dir)}")
            return

        logging.info(f"Filtering {name} (original shape: {adata.shape})...")
        sc.pp.filter_cells(adata, min_genes=200)
        logging.info(f"New shape post-filtering: {adata.shape}")

        adata.obs['dataset'] = name
        adata.obs_names = [f"{name}_{barcode}" for barcode in adata.obs_names]
        adatas.append(adata)

    logging.info("Concatenating datasets...")
    combined_adata = adatas[0].concatenate(adatas[1:], join='outer', index_unique=None)
    
    raw_concat_path = os.path.join(output_dir, 'Concat_dataset_raw.h5ad')
    logging.info(f"Saving concatenated raw dataset to {raw_concat_path}...")
    combined_adata.write(raw_concat_path)

    # --------------------------------------------------------------------------
    # Quality Control and Preprocessing
    # --------------------------------------------------------------------------
    logging.info("Starting preprocessing...")
    
    combined_adata.var_names_make_unique()
    
    logging.info("Calculating QC metrics...")
    combined_adata.var["mito"] = combined_adata.var_names.str.contains("^MT-")
    combined_adata.var["ribo"] = combined_adata.var_names.str.contains("^RP[LS]")
    
    if hasattr(combined_adata.X, 'toarray'):
        combined_adata.var["total_counts"] = np.array(combined_adata.X.sum(0)).flatten()
        combined_adata.var["n_cells"] = np.array((combined_adata.X > 0).sum(0)).flatten()
        combined_adata.var['mean_expr'] = np.array(combined_adata.X.mean(0)).flatten()
        
        combined_adata.obs["total_counts"] = np.array(combined_adata.X.sum(1)).flatten()
        combined_adata.obs["n_genes"] = np.array((combined_adata.X > 0).sum(1)).flatten()
        
        mito_genes = combined_adata.var["mito"].values
        combined_adata.obs["pct_counts_mito"] = (
            np.array(combined_adata.X[:, mito_genes].sum(1)).flatten() / combined_adata.obs["total_counts"] * 100
        )
    else:
        combined_adata.var["total_counts"] = combined_adata.X.sum(0).flatten()
        combined_adata.var["n_cells"] = (combined_adata.X > 0).sum(0).flatten()
        combined_adata.var['mean_expr'] = combined_adata.X.mean(0).flatten()
        
        combined_adata.obs["total_counts"] = combined_adata.X.sum(1).flatten()
        combined_adata.obs["n_genes"] = (combined_adata.X > 0).sum(1).flatten()
        
        mito_genes = combined_adata.var["mito"].values
        combined_adata.obs["pct_counts_mito"] = (
            combined_adata.X[:, mito_genes].sum(1).flatten() / combined_adata.obs["total_counts"] * 100
        )

    combined_adata.obs["log1p_total_counts"] = np.log1p(combined_adata.obs["total_counts"])

    # --------------------------------------------------------------------------
    # Filtering and Normalization
    # --------------------------------------------------------------------------
    logging.info("Filtering based on mitochondrial threshold (<20%)...")
    combined_adata = combined_adata[combined_adata.obs["pct_counts_mito"] < 20].copy()

    logging.info("Normalizing data...")
    combined_adata.layers["counts"] = combined_adata.X.copy()
    combined_adata.layers["lognorm"] = combined_adata.X.copy()
    
    sc.pp.normalize_total(combined_adata, target_sum=1e4, layer="lognorm")
    sc.pp.log1p(combined_adata, layer="lognorm")
    
    if "log1p" in combined_adata.uns:
        del combined_adata.uns["log1p"]

    # --------------------------------------------------------------------------
    # Save Final Processed Object
    # --------------------------------------------------------------------------
    output_path = os.path.join(output_dir, 'Concat_dataset.h5ad')
    logging.info(f"Saving concatenated, QC-filtered dataset to {output_path}...")
    combined_adata.write(output_path)
    
    logging.info("Pipeline Step Complete.")

if __name__ == "__main__":
    main()
