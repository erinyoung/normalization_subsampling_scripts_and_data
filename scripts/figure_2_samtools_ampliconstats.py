import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable

# --- 1. Helper Function to Prepare Heatmap Data ---
def prepare_heatmap_data(filepath):
    df = pd.read_csv(filepath)
    amplicon_cols = df.columns[df.columns.get_loc("Depth") + 1 :]

    df_grouped = (
        df.groupby(["Method", "Depth"])[amplicon_cols].mean().reset_index()
    )
    df_grouped["Depth"] = df_grouped["Depth"].replace("all", np.inf)
    df_grouped["Depth"] = pd.to_numeric(df_grouped["Depth"])
    method_order = {"bbnorm": 0, "seqkit": 1, "raw": 2}
    df_grouped["method_order"] = df_grouped["Method"].map(method_order)
    df_grouped = df_grouped.sort_values(["method_order", "Depth"])

    df_grouped["Method_Depth"] = (
        df_grouped["Method"]
        + "_"
        + df_grouped["Depth"]
        .astype(str)
        .str.replace(".0", "", regex=False)
        .str.replace("inf", "all")
    )

    return df_grouped.set_index("Method_Depth")[amplicon_cols]


# --- 2. Create Figure Grid ---
fig, axes = plt.subplots(2, 2, figsize=(16, 10), constrained_layout=True)
fig.suptitle(
    "Analysis of Amplicon Performance Across Data Reduction Methods",
    fontsize=18,
)

# --- 3. Plot Each Panel ---

# Panel A
heatmap_data_depth = prepare_heatmap_data("data/ampliconstats_fdepth.csv")
heatmap_data_capped = heatmap_data_depth.clip(upper=15000)
sns.heatmap(
    heatmap_data_capped,
    cmap="viridis",
    norm=LogNorm(),
    ax=axes[0, 0],
    cbar_kws={"label": "Mean Read Depth (Log Scale)"},
    rasterized=True,
)
axes[0, 0].set_title(
    "A) Mean Amplicon Read Depth (Log Scale)", loc="left", fontsize=13
)

# Panel B: Pin colorbar explicitly to Panel B's heatmap box
heatmap_data_cov = prepare_heatmap_data("data/ampliconstats_fpcov.csv")
dividerB = make_axes_locatable(axes[0, 1])
caxB = dividerB.append_axes("right", size="3%", pad=0.1)

sns.heatmap(
    heatmap_data_cov,
    cmap="viridis",
    ax=axes[0, 1],
    cbar_ax=caxB,
    cbar_kws={"label": "Percent Coverage"},
    rasterized=True,
)
axes[0, 1].set_title(
    "B) Mean Amplicon Percent Coverage", loc="left", fontsize=13
)

# Panel C
heatmap_data_perc = prepare_heatmap_data("data/ampliconstats_frperc.csv")
sns.heatmap(
    heatmap_data_perc,
    cmap="viridis",
    ax=axes[1, 0],
    cbar_kws={"label": "Percent of Total Reads"},
    rasterized=True,
)
axes[1, 0].set_title(
    "C) Mean Read Percentage per Amplicon", loc="left", fontsize=13
)

for ax in [axes[0, 0], axes[0, 1], axes[1, 0]]:
    ax.set_xlabel("Amplicons (ordered by genome position)", fontsize=10)
    ax.set_xticks([])
    ax.set_ylabel("Method and Depth", fontsize=10)

# Panel D
df_bar = pd.read_csv("data/ampliconstats_fdepth.csv")
amplicon_cols = df_bar.columns[df_bar.columns.get_loc("Depth") + 1 :]
df_bar["Passed (>50X)"] = (df_bar[amplicon_cols] > 50).sum(axis=1)
df_bar["Intermediate (10-50X)"] = (
    (df_bar[amplicon_cols] >= 10) & (df_bar[amplicon_cols] <= 50)
).sum(axis=1)
df_bar["Failed (<10X)"] = (df_bar[amplicon_cols] < 10).sum(axis=1)

agg_df = (
    df_bar.groupby(["Method", "Depth"])[
        ["Passed (>50X)", "Intermediate (10-50X)", "Failed (<10X)"]
    ]
    .mean()
    .reset_index()
)
agg_df["Depth"] = agg_df["Depth"].replace("all", 99999)
agg_df["Depth"] = pd.to_numeric(agg_df["Depth"])
agg_df = agg_df.sort_values(by=["Depth", "Method"])
agg_df["label"] = (
    agg_df["Method"]
    + "_"
    + agg_df["Depth"]
    .astype(str)
    .str.replace(".0", "", regex=False)
    .str.replace("99999", "all")
)
agg_df.set_index("label", inplace=True)

plot_cols = ["Failed (<10X)", "Intermediate (10-50X)", "Passed (>50X)"]
plot_colors = ["black", "gray", "white"]

agg_df[plot_cols].plot(
    kind="bar",
    stacked=True,
    width=0.85,
    ax=axes[1, 1],
    edgecolor="black",
    linewidth=0.8,
    color=plot_colors,
)
axes[1, 1].set_title(
    "D) Average Amplicon Pass/Fail Counts", loc="left", fontsize=13
)
axes[1, 1].set_ylabel("Average Number of Amplicons", fontsize=10)
axes[1, 1].set_xlabel("Method and Depth", fontsize=10)
axes[1, 1].legend(
    title="Depth Category",
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    fontsize=9,
)
plt.setp(
    axes[1, 1].get_xticklabels(),
    rotation=45,
    ha="right",
    rotation_mode="anchor",
)

plt.savefig(
    "Figure2.tiff",
    dpi=300,
    format="tiff",
    pil_kwargs={"compression": "tiff_lzw"},
)