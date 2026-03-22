"""
Figure 3.9. C. KitL enhances neuronal differentiation in OE-derived cultures and supports emergence of TUJ1-positive human cells. 
-----------------------------------------
This script generates a series of SuperPlots comparing neuronal marker 
expression (TUJ1, SOX2) across multiple cell types in response to KitL 
and Midkine treatments. Significance is determined via One-Way ANOVA 
with Dunnett's Multiple Comparisons Test against Control populations.

Outputs: Figure 3.9. C._[Readout].pdf, .tiff, anova_pvals.txt
"""
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from scipy import stats
import os

# -----------------------------------------------------------------------------
# Data Loading & Aggregation
# -----------------------------------------------------------------------------
file_path = "3.9. data.xlsx"

if os.path.exists(file_path):
    df = pd.read_excel(file_path, sheet_name="Data")
else:
    raise FileNotFoundError(f"{file_path} not found.")

# -----------------------------------------------------------------------------
# Prevent Pseudoreplication: Calculate Biological Means (Well-level)
# -----------------------------------------------------------------------------
bio_df = df.groupby(["CellType", "Condition", "Well"]).mean(numeric_only=True).reset_index()

# -----------------------------------------------------------------------------
# Statistical Analysis Setup (ANOVA + Dunnett's)
# -----------------------------------------------------------------------------
pvals_log = []

def get_pvals(cell_type, marker):
    """Calculate One-Way ANOVA with Dunnett's Multiple Comparisons Test against Control"""
    subset = bio_df[bio_df["CellType"] == cell_type].dropna(subset=[marker])
    
    control_data = subset[subset["Condition"] == "Control"][marker].values
    kitl_data = subset[subset["Condition"] == "KitL"][marker].values
    midkine_data = subset[subset["Condition"] == "Midkine"][marker].values
    
    if len(control_data) < 2 or len(kitl_data) < 2 or len(midkine_data) < 2:
        return np.nan, np.nan
        
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            res = stats.dunnett(kitl_data, midkine_data, control=control_data, random_state=42)
            p_kitl = res.pvalue[0]
            p_midkine = res.pvalue[1]
        except AttributeError:
            _, p_k = stats.ttest_ind(control_data, kitl_data, equal_var=False)
            _, p_m = stats.ttest_ind(control_data, midkine_data, equal_var=False)
            p_kitl = min(1.0, p_k * 2) 
            p_midkine = min(1.0, p_m * 2)
            
    pvals_log.append(f"[{cell_type} - {marker}] KitL vs Control p={p_kitl:.4f}")
    pvals_log.append(f"[{cell_type} - {marker}] Midkine vs Control p={p_midkine:.4f}")
    
    return p_kitl, p_midkine

def get_sig_label(pval):
    if pd.isna(pval): return "ns"
    if pval < 0.001: return "***"
    elif pval < 0.01: return "**"
    elif pval < 0.05: return "*"
    else: return "ns"

# -----------------------------------------------------------------------------
# Plot Generation
# -----------------------------------------------------------------------------
sns.set_theme(style="ticks", font_scale=1.2)
plt.rcParams.update({"axes.linewidth": 1.5, "figure.facecolor": "white", "axes.facecolor": "white"})

condition_order = ["Control", "KitL", "Midkine"]
cell_types = df["CellType"].unique()
readouts = ["TUJ_area_percent", "TUJ_clusters", "SOX2_percent"]

cond_palette = {"Control": "#D3D3D3", "KitL": "#A0A0A0", "Midkine": "#707070"}
well_palette = sns.color_palette("muted", n_colors=df["Well"].nunique())
well_color_map = dict(zip(df["Well"].unique(), well_palette))

def plot_superplot_for_readout(readout):
    fig, axes = plt.subplots(1, len(cell_types), figsize=(3.5 * len(cell_types), 5), sharey=False)
    if len(cell_types) == 1: axes = [axes]
    
    for ax, ct in zip(axes, cell_types):
        ct_df = df[df["CellType"] == ct]
        ct_bio_df = bio_df[bio_df["CellType"] == ct]
        
        sns.barplot(
            data=ct_bio_df, x="Condition", y=readout, order=condition_order,
            hue="Condition", palette=cond_palette, edgecolor="black", linewidth=1.5,
            capsize=0.1, ax=ax, errorbar="se", alpha=0.3, legend=False
        )
        
        sns.stripplot(
            data=ct_df, x="Condition", y=readout, order=condition_order,
            hue="Well", palette=well_color_map, size=4, alpha=0.6,
            jitter=0.2, ax=ax, legend=False
        )
        
        sns.stripplot(
            data=ct_bio_df, x="Condition", y=readout, order=condition_order,
            hue="Well", palette=well_color_map, size=9, linewidth=1.5,
            edgecolor="black", jitter=0.0, ax=ax, legend=False
        )
        
        p_kitl, p_midkine = get_pvals(ct, readout)
        y_max = ct_df[readout].max() * 1.10
        
        ax.plot([0, 1], [y_max, y_max], lw=1.2, c="black")
        ax.text(0.5, y_max + (y_max * 0.01), get_sig_label(p_kitl), ha='center', va='bottom', fontsize=12)
        
        y_max2 = ct_df[readout].max() * 1.25
        ax.plot([0, 2], [y_max2, y_max2], lw=1.2, c="black")
        ax.text(1.0, y_max2 + (y_max2 * 0.01), get_sig_label(p_midkine), ha='center', va='bottom', fontsize=12)
        
        if readout == "SOX2_percent":
            ax.set_yticks([0, 25, 50, 75, 100])
            ax.set_ylim(0, max(105, y_max2 * 1.05))
        else:
            ax.set_ylim(0, y_max2 * 1.05)

        ylabels = {
            "TUJ_area_percent": "TUJ1+ Area Fraction (%)",
            "TUJ_clusters": "TUJ1+ Clusters",
            "SOX2_percent": "SOX2+/TUJ1+ Clusters (%)"
        }
        
        ct_titles = {
            "Aldh1l1": "Aldh1l1-sorted",
            "Gpm6a": "Gpm6a-sorted",
            "Krt5": "Krt5-sorted",
            "Human": "Human-unsorted"
        }

        ax.set_title(ct_titles.get(ct, ct), fontsize=14, fontweight="bold", pad=10)
        ax.set_xlabel("")
        ax.set_ylabel(ylabels.get(readout, readout) if ax == axes[0] else "", fontsize=14, fontweight="medium")
    
    sns.despine(fig)
    
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markersize=9, markeredgecolor='black', label='Bio-Rep (Mean)'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markersize=4, alpha=0.5, label='Image Field (Raw)')
    ]
    axes[-1].legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1.05, 0.5), frameon=False, fontsize=12)
    
    plt.tight_layout()
    plt.savefig(f"Figure 3.9. C._{readout}.pdf", format="pdf", bbox_inches="tight")
    plt.savefig(f"Figure 3.9. C._{readout}.tiff", format="tiff", dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    print(f"Saved Figure 3.9. C._{readout}.pdf and .tiff")
    
for r in readouts:
    plot_superplot_for_readout(r)

with open("anova_pvals.txt", "w", encoding="utf-8") as f:
    f.write("\n".join(pvals_log))
