import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# --- 1. Load and Prepare Data ---
# Load the data from CSV file
df = pd.read_csv("data/samtools_coverage.csv")

# Define the order for the x-axis
order = ["10", "30", "50", "100", "200", "300", "400", "500", "1000", "5000", "10000", "all"]
palette = {"bbnorm": "#D55E00", "raw": "#0072B2", "all": "#0072B2", "seqkit": "#009E73"}

# --- 2. Create the Figure Grid ---
# Create a 2x2 grid of plots. `figsize` controls the overall size.
fig, axes = plt.subplots(2, 2, figsize=(15, 12))

# --- 3. Plot Each Panel with Corrected Column Names ---

# Panel A: Covered Bases
sns.boxplot(data=df, 
            x="depth", 
            y="covbases", 
            hue="method",
            order=order, 
            palette=palette, 
            ax=axes[0, 0],
            flierprops={"marker": "o", "markersize": 5, "alpha": 0.3}
            )
axes[0, 0].set_title("A) Covered Bases", fontsize=16, loc='left')
# Remove x-label from top plots
axes[0, 0].set_xlabel("")
axes[0, 0].set_ylabel("Number of Bases Covered")
axes[0, 0].grid(True, linestyle='--', alpha=0.6)

# Panel B: Coverage (%)
sns.boxplot(data=df,
            x="depth", 
            y="coverage", 
            hue="method",
            order=order, 
            palette=palette, 
            ax=axes[0, 1],
            flierprops={"marker": "o", "markersize": 5, "alpha": 0.3}
            )
axes[0, 1].set_title("B) Genome Coverage (%)", fontsize=16, loc='left')
# Remove x-label from top plots
axes[0, 1].set_xlabel("")
axes[0, 1].set_ylabel("Percent of Genome Covered")
# add a reference line for context
axes[0, 1].axhline(95, ls='--', color='red', alpha=0.8)
axes[0, 1].grid(True, linestyle='--', alpha=0.6)

# Panel C: Mean Base Quality
sns.boxplot(data=df, 
            x="depth", 
            y="meanbaseq", 
            hue="method",
            order=order, 
            palette=palette, 
            ax=axes[1, 0],
            flierprops={"marker": "o", "markersize": 5, "alpha": 0.3}
           )
axes[1, 0].set_title("C) Mean Base Quality", fontsize=16, loc='left')
axes[1, 0].set_xlabel("Subsampling Depth Target")
axes[1, 0].set_ylabel("Mean Phred Quality Score")
axes[1, 0].grid(True, linestyle='--', alpha=0.6)

# Panel D: Mean Mapping Quality
sns.boxplot(data=df, 
            x="depth", 
            y="meanmapq", 
            hue="method",
            order=order, 
            palette=palette, 
            ax=axes[1, 1],
            flierprops={"marker": "o", "markersize": 5, "alpha": 0.3}
            )
axes[1, 1].set_title("D) Mean Mapping Quality", fontsize=16, loc='left')
axes[1, 1].set_xlabel("Subsampling Depth Target")
axes[1, 1].set_ylabel("Mean MAPQ Score")
axes[1, 1].grid(True, linestyle='--', alpha=0.6)

# --- 4. Finalize and Save the Figure ---

# Get handles and labels from one of the plots to create a single legend
handles, labels = axes[0, 0].get_legend_handles_labels()

# Hide all the default legends that seaborn automatically creates
for ax in axes.flat:
    ax.get_legend().remove()

# Add a new, single legend and anchor it to the top-right of Panel B
axes[0, 1].legend(handles, labels, title="Method", loc="upper left",
                  bbox_to_anchor=(1, 1))

# Use tight_layout to fix the main plot spacing
plt.tight_layout()

# Manually adjust the figure's right boundary to make space for the legend
fig.subplots_adjust(right=0.87)

# Save the combined figure
plt.savefig("Figure1.png", dpi=300)


