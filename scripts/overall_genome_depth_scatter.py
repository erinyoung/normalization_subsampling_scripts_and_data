import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("data/samtools_coverage.csv")

# Define a color-blind friendly color palette
color_map = {
    'all': '#0173B2',
    'bbnorm': '#DE8F05',
    'seqkit': '#029E73'
}

marker_map = {
    'all': 'o',    # Circle
    'bbnorm': '^', # Triangle
    'seqkit': 's'  # Square
}

# Create a scatter plot
plt.figure(figsize=(10, 6))

# Loop through each method and plot with the specified color and alpha
for method, color in color_map.items():
    method_df = df[df['method'] == method]
    plt.scatter(
        method_df['numreads'], 
        method_df['meandepth'], 
        label=method,
        marker=marker_map[method],  
        color=color, 
        alpha=0.7)

# Add labels and title
plt.xlabel('Number of Reads')
plt.ylabel('Mean Depth')
plt.title('Relationship Between Read Count and Mean Sequencing Depth')
plt.legend()
plt.grid(True)
plt.tight_layout()

# Save the plot to a file
plt.savefig('scatter_meandepth.png', dpi=300)