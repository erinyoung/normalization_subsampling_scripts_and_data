import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# --- 1. Load and Prepare Data ---
df = pd.read_csv("data/pango_lineages.csv")

# --- Prep for Boxplots ---
palette = {'raw': '#0072B2', 'bbnorm': '#D55E00', 'seqkit': '#029E73'}
depth_order = ['10', '30', '50', '100', '200', '300', '400', '500', '1000', '5000', '10000', 'all']
df["ambiguity_content"] = df["qc_notes"].str.extract(r"Ambiguous_content:([0-9.]+)")[0].astype(float)
df['depth'] = df['depth'].astype(str)

# --- Prep for Bar Charts ---
control_df = df[(df['method'] == 'raw') & (df['depth'] == 'all')].copy()
control_df = control_df[['sample', 'lineage']].rename(columns={'lineage': 'control_lineage'})
df_merged = pd.merge(df, control_df, on='sample', how='left')

def classify_outcome_detailed(row):
    lineage = str(row['lineage'])
    control_lineage = str(row['control_lineage'])
    if lineage == control_lineage: return 'Match'
    elif lineage == "Unassigned": return 'Classification Failure'
    else: return 'Mismatch'

comparison_df = df_merged[~((df_merged['method'] == 'raw') & (df_merged['depth'] == 'all'))].copy()
comparison_df['outcome'] = comparison_df.apply(classify_outcome_detailed, axis=1)

summary_table = pd.crosstab(
    index=[comparison_df['method'], comparison_df['depth']],
    columns=comparison_df['outcome']
)
for col in ['Match', 'Mismatch', 'Classification Failure']:
    if col not in summary_table.columns: summary_table[col] = 0
summary_table = summary_table.reset_index()

# --- 2. Create the Figure Grid (2x2) ---
fig, axes = plt.subplots(2, 2, figsize=(20, 16), sharex='col')
fig.suptitle('Impact of Data Reduction on Consensus Quality and Lineage Assignment', fontsize=20, y=0.98)

# --- 3. Plot Each Panel ---

# Panel A: Ambiguity Content (top-left)
sns.boxplot(
    data=df, 
    x='depth', 
    y='ambiguity_content', 
    hue='method', 
    order=depth_order, 
    palette=palette, 
    ax=axes[0, 0],
    flierprops={"marker": "o", "markersize": 5, "alpha": 0.3}
    )
axes[0, 0].set_title('A) Ambiguity Content in Consensus Sequence', loc='left', fontsize=16)
axes[0, 0].set_ylabel('Ambiguity Content (% "N" bases)')
axes[0, 0].grid(True, linestyle='--', alpha=0.6)

# Panel B: Scorpio Support (top-right)
sns.boxplot(
    data=df, 
    x='depth', 
    y='scorpio_support', 
    hue='method', 
    order=depth_order, 
    palette=palette, 
    ax=axes[0, 1],
    flierprops={"marker": "o", "markersize": 5, "alpha": 0.3}
    )
axes[0, 1].set_title('B) Pangolin Scorpio Support', loc='left', fontsize=16)
axes[0, 1].set_ylabel('Scorpio Support Score')
axes[0, 1].grid(True, linestyle='--', alpha=0.6)

# Panel C & D: Lineage Classification Concordance
summary_table['depth'] = pd.to_numeric(summary_table['depth'])
summary_table = summary_table.sort_values(by=['depth'])
plot_df = summary_table.set_index(['method', 'depth'])
plot_cols = ['Classification Failure', 'Mismatch', 'Match']

# Plot BBNorm on bottom-left
bbnorm_df = plot_df.loc['bbnorm']
bbnorm_df[plot_cols].plot(
    kind='bar', 
    stacked=True, 
    ax=axes[1, 0], 
    color=['black', 'gray', 'white'],
    edgecolor='black', 
    legend=False, width=0.8
    )
bars = axes[1, 0].patches
axes[1, 0].set_title('C) BBNorm Lineage Outcomes', loc='left', fontsize=16)
axes[1, 0].set_ylabel('Number of Samples')
axes[1, 0].set_xlabel('Subsampling Depth Target')

# Plot SeqKit on bottom-right
seq_df = plot_df.loc['seqkit']
seq_df[plot_cols].plot(kind='bar', 
                       stacked=True, 
                       ax=axes[1, 1], 
                       color=['black', 'grey', 'white'], 
                       edgecolor='black',
                       legend=False, 
                       width=0.8
                       )
bars = axes[1, 1].patches
axes[1, 1].set_title('D) SeqKit Lineage Outcomes', loc='left', fontsize=16)
axes[1, 1].set_xlabel('Subsampling Depth Target')
axes[1, 1].set_ylabel('')

# --- 4. Finalize Legends and Layout ---
# Remove the default legends from the boxplots
axes[0, 0].get_legend().remove()
handles, labels = axes[0, 1].get_legend_handles_labels()
axes[0, 1].get_legend().remove()

# Add a shared legend for the boxplots, anchored to the outside of Panel B
axes[0, 1].legend(handles, labels, title='Method', loc='upper left', bbox_to_anchor=(1.02, 1))

# Add a shared legend for the bar charts
legend_elements = [
    Patch(facecolor='black', edgecolor='black', label='Failure'),
    Patch(facecolor='gray',  edgecolor='black', label='Mismatch'),
    Patch(facecolor='white', edgecolor='black', label='Match')
]
axes[1, 1].legend(handles=legend_elements, title='Outcome', loc='upper left', bbox_to_anchor=(1.02, 1))

# Adjust layout to make space for the legends
plt.tight_layout(rect=[0, 0, 0.9, 0.96])
plt.savefig("Figure3.png", dpi=300)

plt.show()