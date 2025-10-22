import pandas as pd
import io
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.read_csv("data/pango_lineages.csv")

boxplot_color_map = {
    'raw': '#0173B2',    # Blue
    'bbnorm': '#DE8F05', # Orange
    'seqkit': '#029E73'  # Green
}

# The 'depth' column has mixed types, convert to string for filtering
df['depth_str'] = df['depth'].astype(str)

df["ambiguity_content"] = (
    df["qc_notes"]
    .str.extract(r"Ambiguous_content:([0-9.]+)")[0]
    .astype(float)
)

# Filter for the specific depth levels of interest
desired_depths = ['10', '30', '50','100', '200', '300', '400', '500', '1000', '5000', '10000', 'all']
plot_df = df[df['depth_str'].isin(desired_depths)].copy()

# Set a specific order for the plot categories
depth_order = ['10', '30', '50','100', '200', '300', '400', '500', '1000', '5000', '10000', 'all']

# --- Plot ---
plt.figure(figsize=(12, 7))
sns.boxplot(
    data=plot_df, 
    x='depth_str', 
    y='ambiguity_content', 
    hue='method', 
    order=depth_order,
    palette=boxplot_color_map
    )
plt.title('Ambiguity of Content by Method and Subsampling Depth', fontsize=16)
plt.xlabel('Subsampling Depth Target', fontsize=12)
plt.ylabel('Ambiguity Content', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.6)
plt.savefig('boxplot_ambiguity_content.png', dpi=300)
plt.savefig('boxplot_ambiguity_content.svg', dpi=300)
plt.close()