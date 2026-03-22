"""
Figure 3.7. C. Aldh1l1-enriched cultures exhibit greater neurogenic capacity than Krt5+ basal cultures
-----------------------------------------
This script generates a SuperPlot comparing the area fraction of neuronal 
markers (TUJ1, TUBB4) between Aldh1l1-GFP and Krt5-TdT cell pools.
Significance is determined via a two-sided Welch's t-test on biological means.

Outputs: Figure 3.7 C.pdf, Figure 3.7 C.tiff
"""
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from scipy import stats

# -----------------------------------------------------------------------------
# Data Definition & Processing
# -----------------------------------------------------------------------------
data_records = [
    # Aldh1l1
    {"Group": "Aldh1l1-GFP", "Replicate": "Rep1", "Field": "Field1", "Tuj1": 12.7, "B4": 15.2},
    {"Group": "Aldh1l1-GFP", "Replicate": "Rep1", "Field": "Field2", "Tuj1": 6.5,  "B4": 17.2},
    {"Group": "Aldh1l1-GFP", "Replicate": "Rep1", "Field": "Field3", "Tuj1": 7.04, "B4": 13.1},
    {"Group": "Aldh1l1-GFP", "Replicate": "Rep2", "Field": "Field1", "Tuj1": 5.32, "B4": 20.8},
    {"Group": "Aldh1l1-GFP", "Replicate": "Rep2", "Field": "Field2", "Tuj1": 6.25, "B4": 20.3},
    {"Group": "Aldh1l1-GFP", "Replicate": "Rep2", "Field": "Field3", "Tuj1": 6.21, "B4": 18.4},
    # KT
    {"Group": "Krt5-TdT", "Replicate": "Rep1", "Field": "Field1", "Tuj1": 0.0,  "B4": 24.7},
    {"Group": "Krt5-TdT", "Replicate": "Rep1", "Field": "Field2", "Tuj1": 0.0,  "B4": 27.3},
    {"Group": "Krt5-TdT", "Replicate": "Rep1", "Field": "Field3", "Tuj1": 0.0,  "B4": 25.5},
    {"Group": "Krt5-TdT", "Replicate": "Rep2", "Field": "Field1", "Tuj1": 0.01, "B4": 25.9},
    {"Group": "Krt5-TdT", "Replicate": "Rep2", "Field": "Field2", "Tuj1": 0.0,  "B4": 36.6},
    {"Group": "Krt5-TdT", "Replicate": "Rep2", "Field": "Field3", "Tuj1": 0.0,  "B4": 28.9},
]

df = pd.DataFrame(data_records)

# -----------------------------------------------------------------------------
# SuperPlot Aggregation: Calculate Biological Means
# -----------------------------------------------------------------------------
bio_df = df.groupby(["Group", "Replicate"]).mean(numeric_only=True).reset_index()

# -----------------------------------------------------------------------------
# Statistical Analysis (Welch's t-test on Biological Means)
# -----------------------------------------------------------------------------
def run_stats(marker):
    group1 = bio_df[bio_df["Group"] == "Aldh1l1-GFP"][marker]
    group2 = bio_df[bio_df["Group"] == "Krt5-TdT"][marker]
    t_stat, p_val = stats.ttest_ind(group1, group2, equal_var=False)
    return p_val

pval_tuj1 = run_stats("Tuj1")
pval_b4 = run_stats("B4")

def get_sig_label(pval):
    if np.isnan(pval): return "ns"
    if pval < 0.001: return "***"
    elif pval < 0.01: return "**"
    elif pval < 0.05: return "*"
    else: return "ns"

sig_tuj1 = get_sig_label(pval_tuj1)
sig_b4 = get_sig_label(pval_b4)

# -----------------------------------------------------------------------------
# Plot Styling
# -----------------------------------------------------------------------------
sns.set_theme(style="ticks", font_scale=1.2)
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.linewidth": 1.5,
    "axes.edgecolor": "black",
    "axes.labelcolor": "black",
    "xtick.color": "black",
    "ytick.color": "black",
    "xtick.major.width": 1.5,
    "ytick.major.width": 1.5,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
})

# -----------------------------------------------------------------------------
# Plot Generation
# -----------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(8.5, 5), sharey=True)

rep_palette = {"Rep1": "#5C82A6", "Rep2": "#D38D5F"}

def draw_superpanel(ax, marker, title, pval, sig_label):
    sns.barplot(
        data=bio_df,
        x="Group",
        y=marker,
        color="#E0E0E0",
        edgecolor="black",
        linewidth=1.5,
        capsize=0.1,
        ax=ax,
        errorbar="se",
        alpha=0.5
    )
    
    sns.stripplot(
        data=df,
        x="Group",
        y=marker,
        hue="Replicate",
        palette=rep_palette,
        size=4,
        alpha=0.5,
        jitter=0.15,
        ax=ax,
        legend=False
    )
    
    sns.stripplot(
        data=bio_df,
        x="Group",
        y=marker,
        hue="Replicate",
        palette=rep_palette,
        size=9,
        linewidth=1.5,
        edgecolor="black",
        jitter=0.0,
        ax=ax,
        legend=False
    )
    
    y_max = df[marker].max()
    y_sig = y_max * 1.05
    y_line_offset = y_max * 0.03
    
    ax.plot([0, 1], [y_sig, y_sig], lw=1.5, c="black")
    ax.text(0.5, y_sig + y_line_offset, sig_label, ha='center', va='bottom', color="black", fontsize=14)
    
    ax.set_title(title, fontsize=14, fontweight="bold", pad=10)
    ax.set_xlabel("", fontsize=12)
    ax.set_ylabel("Area Fraction (%)", fontsize=14, fontweight="medium")

draw_superpanel(axes[0], "Tuj1", "TUJ1", pval_tuj1, sig_tuj1)
draw_superpanel(axes[1], "B4", "TUBB4", pval_b4, sig_b4)
axes[1].set_ylabel("")

y_global_max = df[["Tuj1", "B4"]].max().max()
axes[0].set_ylim(0, y_global_max * 1.25)

sns.despine(fig, top=True, right=True)

# -----------------------------------------------------------------------------
# Construct Custom Legend
# -----------------------------------------------------------------------------
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches

legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markersize=9, markeredgecolor='black', label='Bio-Rep (Mean)'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markersize=4, alpha=0.5, label='Image Field (Raw)'),
]
axes[1].legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1.10, 0.5), frameon=False, fontsize=12)

plt.tight_layout()

# -----------------------------------------------------------------------------
# Output & Save
# -----------------------------------------------------------------------------
plt.savefig("Figure 3.7 C.pdf", format="pdf", bbox_inches="tight")
plt.savefig("Figure 3.7 C.tiff", format="tiff", dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})

print("Tuj1 Welch's t-test p-value (N=2):", pval_tuj1)
print("B4 Welch's t-test p-value (N=2):", pval_b4)
