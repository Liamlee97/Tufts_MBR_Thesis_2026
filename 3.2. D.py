"""
Figure 3.2. D. Murine and human epithelial cultures exhibit divergent differentiation outcomes under matched conditions
-----------------------------------------
This script generates publication-quality bar and strip plots comparing the area 
fraction of neuronal markers (TUJ1, TUBB4) between Mouse and Human cultures. 
Significance is determined via a two-sided Welch's t-test.

Outputs: Figure 3.2 D.pdf, Figure 3.2 D.tiff
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
    # Mouse Data
    {"Species": "Mouse", "Sample": "VC OG", "Tuj1": 61.6, "B4": 16.7},
    {"Species": "Mouse", "Sample": "VC361", "Tuj1": 16.1, "B4": 4.39},
    {"Species": "Mouse", "Sample": "VC365", "Tuj1": 30.1, "B4": 13.8},
    {"Species": "Mouse", "Sample": "VC365 #2", "Tuj1": 25.5, "B4": 24.1},
    {"Species": "Mouse", "Sample": "VC411", "Tuj1": 11.859, "B4": 16.6},
    {"Species": "Mouse", "Sample": "LL CD1 w2", "Tuj1": 30.7, "B4": 14.4},
    # Human Data
    {"Species": "Human", "Sample": "HuN8", "Tuj1": 1.61, "B4": 30.5},
    {"Species": "Human", "Sample": "HuN4", "Tuj1": 0.0, "B4": 12.4},
    {"Species": "Human", "Sample": "HuN7", "Tuj1": 0.0, "B4": 8.96},
]

df = pd.DataFrame(data_records)

# -----------------------------------------------------------------------------
# Statistical Analysis (Welch's t-test)
# -----------------------------------------------------------------------------
def run_stats(marker):
    mouse_data = df[df["Species"] == "Mouse"][marker]
    human_data = df[df["Species"] == "Human"][marker]
    # Welch's t-test (assumes unequal variances, preferred for unequal sample sizes)
    t_stat, p_val = stats.ttest_ind(mouse_data, human_data, equal_var=False)
    return p_val

pval_tuj1 = run_stats("Tuj1")
pval_b4 = run_stats("B4")

def get_sig_label(pval):
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
fig, axes = plt.subplots(1, 2, figsize=(6.5, 5), sharey=True)

palette = {"Mouse": "#D3D3D3", "Human": "#A0A0A0"} 
dot_color = "#333333"

def draw_panel(ax, marker, title, pval, sig_label):
    sns.barplot(
        data=df,
        x="Species",
        y=marker,
        hue="Species",
        palette=palette,
        edgecolor="black",
        linewidth=1.5,
        capsize=0.1,
        ax=ax,
        errorbar="se",
        legend=False
    )
    
    sns.stripplot(
        data=df,
        x="Species",
        y=marker,
        color=dot_color,
        size=7,
        jitter=0.15,
        alpha=0.7,
        ax=ax
    )
    
    y_max = df[marker].max()
    y_sig = y_max * 1.1
    y_line_offset = y_max * 0.03
    
    ax.plot([0, 1], [y_sig, y_sig], lw=1.5, c="black")
    ax.text(0.5, y_sig + y_line_offset, sig_label, ha='center', va='bottom', color="black", fontsize=14)
    
    ax.set_title(title, fontsize=14, fontweight="bold", pad=10)
    ax.set_xlabel("", fontsize=12)
    ax.set_ylabel("Area Fraction (%)", fontsize=14, fontweight="medium")

draw_panel(axes[0], "Tuj1", "TUJ1", pval_tuj1, sig_tuj1)
draw_panel(axes[1], "B4", "TUBB4", pval_b4, sig_b4)
axes[1].set_ylabel("")

y_global_max = df[["Tuj1", "B4"]].max().max()
axes[0].set_ylim(0, y_global_max * 1.35)

sns.despine(fig, top=True, right=True)

plt.tight_layout()

# -----------------------------------------------------------------------------
# Output & Save
# -----------------------------------------------------------------------------
plt.savefig("Figure 3.2 D.pdf", format="pdf", bbox_inches="tight")
plt.savefig("Figure 3.2 D.tiff", format="tiff", dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})

print(f"Tuj1: Mouse vs Human Welch's t-test p-value: {pval_tuj1:.4f}")
print(f"B4: Mouse vs Human Welch's t-test p-value: {pval_b4:.4f}")
