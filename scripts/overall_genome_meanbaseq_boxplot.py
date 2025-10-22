import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.read_csv("data/samtools_coverage.csv")

boxplot_color_map = {
    'all': '#0173B2',    # Blue
    'bbnorm': '#DE8F05', # Orange
    'seqkit': '#029E73'  # Green
}

# The 'depth' column has mixed types, convert to string for filtering
df['depth_str'] = df['depth'].astype(str)

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
    y='meanbaseq', 
    hue='method', 
    order=depth_order,
    palette=boxplot_color_map
    )
plt.title('Distribution of Mean Base Quality by Method and Subsampling Depth', fontsize=16)
plt.xlabel('Subsampling Depth Target', fontsize=12)
plt.ylabel('Mean Base Quality', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.6)
plt.savefig('boxplot_meanbaseq.png', dpi=300)
plt.savefig('boxplot_meanbaseq.svg', dpi=300)
plt.close()