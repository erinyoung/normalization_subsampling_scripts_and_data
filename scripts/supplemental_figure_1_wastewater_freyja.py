import pandas as pd
import matplotlib.pyplot as plt
import os

# 1. Load the data
file_path = 'data/wastewater_freyja.tsv'
df = pd.read_csv(file_path, sep='\t')

# 2. Helper function to parse sample names for grouping
def parse_sample_name(raw_name):
    """Parses base sample name, method, and depth from the filename."""
    raw_name = raw_name.replace('_variants.tsv', '')
    if '_bbnorm_' in raw_name:
        base, depth = raw_name.split('_bbnorm_')
        return base, 'BBNorm', int(depth)
    elif '_seqkit_' in raw_name:
        base, depth = raw_name.split('_seqkit_')
        return base, 'SeqKit', int(depth)
    else:
        # Assign infinite depth to 'Original' so it sorts to the very top
        return raw_name, 'Original', float('inf')

# Apply parsing to create new grouping columns
df['BaseSample'], df['Method'], df['Depth'] = zip(*df['Unnamed: 0'].apply(parse_sample_name))

# 3. Generate a separate plot for each base wastewater sample
base_samples = df['BaseSample'].unique()
os.makedirs('supplemental_figures', exist_ok=True)

for base_sample in base_samples:
    # Isolate data for just this one wastewater sample
    sample_df = df[df['BaseSample'] == base_sample].copy()
    
    # Sort subsamples by Depth (high to low) and Method
    sample_df = sample_df.sort_values(by=['Depth', 'Method'], ascending=[False, True])
    
    # Calculate the top 19 most abundant lineages FOR THIS SPECIFIC SAMPLE
    sample_totals = {}
    for index, row in sample_df.iterrows():
        if pd.isna(row['lineages']) or pd.isna(row['abundances']):
            continue
        lins = str(row['lineages']).split()
        abunds = str(row['abundances']).split()
        
        for lin, abund in zip(lins, abunds):
            try:
                sample_totals[lin] = sample_totals.get(lin, 0) + float(abund)
            except ValueError:
                continue
                
    # Sort and keep the top 19 specific to this isolate
    top_19_lins = sorted(sample_totals, key=sample_totals.get, reverse=True)[:19]
    
    parsed_data = []
    y_labels = []
    
    for index, row in sample_df.iterrows():
        if pd.isna(row['lineages']):
            continue
        
        method = row['Method']
        depth = row['Depth']
        
        # Format the Y-axis label
        if method == 'Original':
            label = 'Original'
        else:
            label = f"{int(depth)}X {method}"
            
        y_labels.append(label)
        
        lins = str(row['lineages']).split()
        abunds = str(row['abundances']).split()
        
        sample_dict = {}
        other_val = 0
        
        for lin, abund in zip(lins, abunds):
            try:
                val = float(abund)
                # Categorize into this sample's Top 19 or "Other"
                if lin in top_19_lins:
                    sample_dict[lin] = sample_dict.get(lin, 0) + val
                else:
                    other_val += val
            except ValueError:
                continue
                
        if other_val > 0:
            sample_dict['Other'] = other_val
            
        parsed_data.append(sample_dict)
        
    if not parsed_data:
        continue
        
    # Build plot dataframe
    plot_df = pd.DataFrame(parsed_data, index=y_labels).fillna(0)
    
    # Reorder columns to ensure "Other" is always at the far right
    cols = [c for c in plot_df.columns if c != 'Other']
    if 'Other' in plot_df.columns:
        cols.append('Other')
    plot_df = plot_df[cols]
    
    # Reverse the dataframe so that the highest depths sit at the top of the chart
    plot_df = plot_df.iloc[::-1]
    
    # 4. Plotting Configuration
    fig_height = max(6, len(plot_df) * 0.3)
    fig, ax = plt.subplots(figsize=(10, fig_height))
    
    # Use tab20 to distinctively color separate categories
    plot_df.plot(
        kind='barh', 
        stacked=True, 
        ax=ax, 
        colormap='tab20', 
        edgecolor='white', 
        linewidth=0.3,
        width=0.85
    )
    
    # Axis and Label Formatting
    ax.set_xlabel('Variant Prevalence', fontsize=12, fontweight='bold')
    ax.set_ylabel('Subsampling Method & Depth', fontsize=12, fontweight='bold')
    ax.set_title(f'Lineage Composition: {base_sample}', fontsize=14, fontweight='bold')
    ax.set_xlim(0, 1.0)
    
    ax.tick_params(axis='y', labelsize=9) 
    
    # Format external legend
    ax.legend(
        title='Lineage',
        title_fontsize='11',
        bbox_to_anchor=(1.02, 1), 
        loc='upper left', 
        frameon=False,
        fontsize=9.5
    )
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    output_filename = f'supplemental_figures/wastewater_lineage_{base_sample}.png'
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    plt.close(fig)