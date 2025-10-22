import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

df = pd.read_csv("data/ampliconstats_fpcov.csv")

# Keep only amplicon columns
amplicon_cols = df.columns[5:]

# Group by Method and Depth, take the mean across samples
df_grouped = df.groupby(['Method', 'Depth'])[amplicon_cols].mean().reset_index()

# Create a Method_Depth label for the heatmap rows
df_grouped['Method_Depth'] = df_grouped['Method'] + df_grouped['Depth'].astype(str)

# Sort rows: bbnorm → seqkit → raw (control at bottom)
df_grouped['Depth_numeric'] = pd.to_numeric(df_grouped['Depth'], errors='coerce')
method_order = {'bbnorm': 0, 'seqkit': 1, 'raw': 2}
df_grouped['method_order'] = df_grouped['Method'].map(method_order)
df_grouped = df_grouped.sort_values(['method_order', 'Depth_numeric'])

# Set Method_Depth as index for heatmap
heatmap_data = df_grouped.set_index('Method_Depth')[amplicon_cols]

# Plot the heatmap
plt.figure(figsize=(20, 12))
sns.heatmap(heatmap_data, cmap="viridis", cbar_kws={'label': 'FPCOV'})
plt.title("Mean Ampilcon Percent Coverage by Method and Depth")
plt.xlabel("Amplicons")
plt.ylabel("Method / Depth")
plt.tight_layout()
plt.savefig("samtools_ampliconstats_fpcov_heatmap.png", dpi=300)
plt.savefig("samtools_ampliconstats_fpcov_heatmap.svg", dpi=300)
plt.close()
