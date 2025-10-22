import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

# --- 1. Helper Function to Prepare Heatmap Data ---
# This function handles the repetitive data loading and shaping for the heatmaps.
def prepare_heatmap_data(filepath):
    """Loads a CSV, calculates the mean across samples, and sorts for plotting."""
    df = pd.read_csv(filepath)
    # Correctly identify amplicon columns (they start after the 'Depth' column)
    amplicon_cols = df.columns[df.columns.get_loc('Depth') + 1:]
    
    # Group by Method and Depth, take the mean across samples
    df_grouped = df.groupby(['Method', 'Depth'])[amplicon_cols].mean().reset_index()

    # Create a consistent sorting order
    df_grouped['Depth'] = df_grouped['Depth'].replace('all', np.inf)
    df_grouped['Depth'] = pd.to_numeric(df_grouped['Depth'])
    method_order = {'bbnorm': 0, 'seqkit': 1, 'raw': 2}
    df_grouped['method_order'] = df_grouped['Method'].map(method_order)
    df_grouped = df_grouped.sort_values(['method_order', 'Depth'])

    # Create labels for the y-axis
    df_grouped['Method_Depth'] = df_grouped['Method'] + "_" + df_grouped['Depth'].astype(str).str.replace(".0", "", regex=False).str.replace("inf", "all")
    
    heatmap_data = df_grouped.set_index('Method_Depth')[amplicon_cols]
    return heatmap_data

# --- 2. Create the Figure Grid ---
# Create a 2x2 grid. `constrained_layout=True` helps with spacing.
fig, axes = plt.subplots(2, 2, figsize=(24, 18), constrained_layout=True)
fig.suptitle("Analysis of Amplicon Performance Across Data Reduction Methods", fontsize=24)


# --- 3. Plot Each Panel ---

# Panel A: Mean Amplicon Read Depth (with Logarithmic Scaling)
heatmap_data_depth = prepare_heatmap_data("data/ampliconstats_fdepth.csv")
# Cap the maximum value to improve color scale readability
heatmap_data_capped = heatmap_data_depth.clip(upper=15000)
sns.heatmap(
    heatmap_data_capped,
    cmap="viridis",        # A good, vibrant palette for log scale
    norm=LogNorm(),        # Apply logarithmic scaling
    ax=axes[0, 0],
    cbar_kws={'label': 'Mean Read Depth (Log Scale)'}
)
axes[0, 0].set_title("A) Mean Amplicon Read Depth (Log Scale)", loc='left', fontsize=16)


# Panel B: Mean Amplicon Percent Coverage
heatmap_data_cov = prepare_heatmap_data("data/ampliconstats_fpcov.csv")
sns.heatmap(heatmap_data_cov, cmap="viridis", ax=axes[0, 1], cbar_kws={'label': 'Percent Coverage'})
axes[0, 1].set_title("B) Mean Amplicon Percent Coverage", loc='left', fontsize=16)


# Panel C: Mean Read Percentage Distribution
heatmap_data_perc = prepare_heatmap_data("data/ampliconstats_frperc.csv")
sns.heatmap(heatmap_data_perc, cmap="viridis", ax=axes[1, 0], cbar_kws={'label': 'Percent of Total Reads'})
axes[1, 0].set_title("C) Mean Read Percentage per Amplicon", loc='left', fontsize=16)

# Simplify x-axis labels for all heatmaps
for ax in [axes[0, 0], axes[0, 1], axes[1, 0]]:
    ax.set_xlabel("Amplicons (ordered by genome position)")
    ax.set_xticks([]) # Hide the long, unreadable labels
    ax.set_ylabel("Method and Depth")


# Panel D: Amplicon Pass/Fail Bar Chart (with intuitive colors)
df_bar = pd.read_csv("data/ampliconstats_fdepth.csv")
amplicon_cols = df_bar.columns[df_bar.columns.get_loc('Depth') + 1:]
df_bar['Passed (>50X)'] = (df_bar[amplicon_cols] > 50).sum(axis=1)
df_bar['Intermediate (10-50X)'] = ((df_bar[amplicon_cols] >= 10) & (df_bar[amplicon_cols] <= 50)).sum(axis=1)
df_bar['Failed (<10X)'] = (df_bar[amplicon_cols] < 10).sum(axis=1)
agg_df = df_bar.groupby(['Method', 'Depth'])[['Passed (>50X)', 'Intermediate (10-50X)', 'Failed (<10X)']].mean().reset_index()

# Reorder for grouped comparison
agg_df['Depth'] = agg_df['Depth'].replace('all', 99999) 
agg_df['Depth'] = pd.to_numeric(agg_df['Depth'])
agg_df = agg_df.sort_values(by=['Depth', 'Method'])
agg_df['label'] = agg_df['Method'] + "_" + agg_df['Depth'].astype(str).str.replace(".0", "", regex=False).str.replace("99999", "all")
agg_df.set_index('label', inplace=True)

# Define column order for plotting (puts Failed at the bottom) and matching colors
plot_cols = ['Failed (<10X)', 'Intermediate (10-50X)', 'Passed (>50X)']
plot_colors = ['black', 'gray', 'white']

# Plot the bar chart on the specified subplot
agg_df[plot_cols].plot(
    kind='bar',
    stacked=True,
    ax=axes[1, 1],
    edgecolor='black',
    linewidth=1,
    color=plot_colors
)
axes[1, 1].set_title('D) Average Amplicon Pass/Fail Counts', loc='left', fontsize=16)
axes[1, 1].set_ylabel('Average Number of Amplicons')
axes[1, 1].set_xlabel('Method and Depth')
axes[1, 1].legend(title='Depth Category', bbox_to_anchor=(1.02, 1), loc='upper left')
plt.setp(axes[1, 1].get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")


# --- 4. Finalize and Save the Figure ---
plt.savefig("Figure2", dpi=300, bbox_inches='tight')

plt.show()