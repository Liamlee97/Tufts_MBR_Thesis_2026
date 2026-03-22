"""
Figure 3.6. C. GPM6A-based sorting enriches OE-associated cells but does not isolate a pure basal population.
-----------------------------------------
This script generates publication-quality bar and strip plots quantifying the 
percentage of KRT14-positive cells from DAPI-detected nuclei per well.

Outputs: Figure 3.6 C.pdf, Figure 3.6 C.tiff
"""
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# -----------------------------------------------------------------------------
# Data Definition & Processing
# -----------------------------------------------------------------------------
# Data represents DAPI-detected nuclei per well and Krt14-positive subset
data_records = [
    {
        "Well": "Negative sort",
        "Total_DAPI": 5255,
        "Krt14_Positive": 10
    },
    {
        "Well": "Positive sort",
        "Total_DAPI": 5644,
        "Krt14_Positive": 59
    }
]

df = pd.DataFrame(data_records)

df["Percent_Positive"] = (df["Krt14_Positive"] / df["Total_DAPI"]) * 100

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
fig, ax = plt.subplots(figsize=(4, 5))

bar_color = "#D3D3D3" 
dot_color = "#333333" 

sns.barplot(
    data=df,
    x="Well",
    y="Percent_Positive",
    color=bar_color,
    edgecolor="black",
    linewidth=1.5,
    capsize=0.1,
    ax=ax,
    errorbar=None 
)

sns.stripplot(
    data=df,
    x="Well",
    y="Percent_Positive",
    color=dot_color,
    size=8,
    jitter=False, 
    ax=ax
)

for idx, row in df.iterrows():
    pct_text = f"{row['Percent_Positive']:.2f}%"
    count_text = f"({row['Krt14_Positive']}/{row['Total_DAPI']})"
    
    y_offset = df["Percent_Positive"].max() * 0.05
    
    ax.text(
        idx, 
        row['Percent_Positive'] + y_offset, 
        f"{pct_text}\n{count_text}", 
        ha="center", 
        va="bottom", 
        fontsize=10, 
        color="black"
    )

ax.set_ylabel("% KRT14-positive cells (of DAPI+)", fontsize=14, fontweight="medium")
ax.set_xlabel("", fontsize=12)
ax.set_title("KRT14-positive Cell Frequency", fontsize=14, fontweight="bold", pad=15)

ax.set_ylim(0, df["Percent_Positive"].max() * 1.3)

sns.despine(ax=ax, top=True, right=True)

plt.tight_layout()

# -----------------------------------------------------------------------------
# Output & Save
# -----------------------------------------------------------------------------
plt.savefig("Figure 3.6 C.pdf", format="pdf", bbox_inches="tight")

plt.savefig("Figure 3.6 C.tiff", format="tiff", dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
