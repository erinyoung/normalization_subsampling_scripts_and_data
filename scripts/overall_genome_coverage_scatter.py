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

plt.axvline(x=1980,    color='grey', linestyle='--', alpha = 0.5, label='10X')
plt.axvline(x=19803,   color='grey', linestyle='--', alpha = 0.6, label='100X')
plt.axvline(x=198033,  color='grey', linestyle='--', alpha = 0.8, label='1000X')
plt.axvline(x=1980331, color='grey', linestyle='--', alpha = 1.0, label='10000X')

# Loop through each method and plot with the specified color and alpha
for method, color in color_map.items():
    method_df = df[df['method'] == method]
    plt.scatter(
        method_df['numreads'], 
        method_df['coverage'], 
        label=method, 
        color=color, 
        marker=marker_map[method], 
        alpha=0.7
        )

plt.xscale('log')

# Add labels and title
plt.xlabel('Number of Reads')
plt.ylabel('Coverage')
plt.title('Relationship Between Read Count and Coverage')
plt.legend()
plt.grid(True)
plt.tight_layout()

# Save the plot to a file
plt.savefig('scatter_coverage.png', dpi=300)