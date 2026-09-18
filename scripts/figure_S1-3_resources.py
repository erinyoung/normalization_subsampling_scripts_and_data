import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import math
import string

# --- 1. Data Preparation Functions ---

def parse_col_name(col_name):
    """Parses the complex column names from the wide CSV file."""
    parts = col_name.split('_')
    method, step, depth, resource = "", "", "", ""
    if "bbnorm" in col_name: method = "bbnorm"
    elif "seqkit" in col_name: method = "seqkit"
    elif "raw" in col_name: method = "raw"
    resource = parts[-1]
    if "bbnorm_bbnorm" in col_name:
        step = "bbnorm"
        depth = parts[-2]
    elif "seqkit_seqkit" in col_name:
        step = f"seqkit_{parts[3]}"
        depth = parts[2]
    elif "raw" in col_name:
        step = "_".join(parts[1:-2])
        depth = parts[-2]
    else:
        step = "_".join(parts[1:-2])
        depth = parts[-2]
    num_depth_mapping = {
        '990': '10', '2970': '30', '4951': '50', '9902': '100', '14852': '200',
        '19803': '300', '29705': '400', '49508': '500', '99017': '1000',
        '495083': '5000', '990166': '10000'
    }
    if method == "seqkit" and "seqkit" in step:
        depth = num_depth_mapping.get(depth, depth)
    return method, step, depth, resource

def time_to_seconds(time_str):
    """Converts H:M:S or M:S time strings to total seconds."""
    parts = str(time_str).split(':')
    seconds = 0
    if len(parts) == 3:
        seconds = int(parts[0]) * 3600 + int(parts[1]) * 60 + float(parts[2])
    elif len(parts) == 2:
        seconds = int(parts[0]) * 60 + float(parts[1])
    elif time_str:
        try: seconds = float(time_str)
        except (ValueError, TypeError): seconds = 0
    return seconds

def load_and_reshape_data(filepath):
    """Loads the wide CSV and transforms it into a clean, long-format DataFrame."""
    df = pd.read_csv(filepath)
    df_long = df.melt(id_vars=['sample'], var_name='metric', value_name='value')
    parsed_metrics = df_long['metric'].apply(parse_col_name).apply(pd.Series)
    parsed_metrics.columns = ['method', 'step', 'depth', 'resource']
    df_long = pd.concat([df_long.drop(columns='metric'), parsed_metrics], axis=1)
    is_time = df_long['resource'] == 'time'
    df_long.loc[is_time, 'value'] = df_long.loc[is_time, 'value'].apply(time_to_seconds)
    df_long['value'] = pd.to_numeric(df_long['value'], errors='coerce').fillna(0)
    return df_long

# --- 2. Reusable Multi-Panel Plotting Function ---

def create_multi_panel_supplemental_figure(df, resource_type, output_dir):
    """Generates a single multi-panel figure for a given resource type."""
    
    resource_df = df[df['resource'] == resource_type].copy()
    resource_df['plot_step'] = resource_df['step'].str.replace(r'_R[12]$', '', regex=True)
    
    # *** NEW: Define the desired order for the panels ***
    panel_order = [
        'seqkit', 'bbnorm', 'fastp', 'bbmap', 'sort', 'ivar_trim', 
        'sort2', 'ivar_consensus_mpileup', 'ivar_consensus_consensus'
    ]
    # Filter the unique steps found in the data to match the desired order
    all_steps_in_data = resource_df['plot_step'].unique()
    steps_to_plot = [step for step in panel_order if step in all_steps_in_data]

    # Determine grid size
    num_steps = len(steps_to_plot)
    num_cols = 3
    num_rows = math.ceil(num_steps / num_cols)

    fig, axes = plt.subplots(num_rows, num_cols, figsize=(num_cols * 7, num_rows * 5))
    fig.suptitle(f'Supplemental Figure: Per-Step {resource_type.capitalize()} Usage', fontsize=20, y=1.02)
    
    axes = axes.flatten()
    
    palette = {'raw': '#0072B2', 'bbnorm': '#D55E00', 'seqkit': '#029E73'}
    depth_order = ['10', '30', '50', '100', '200', '300', '400', '500', '1000', '5000', '10000', 'all']
    panel_labels = string.ascii_uppercase

    for i, step in enumerate(steps_to_plot):
        ax = axes[i]
        plot_df = resource_df[resource_df['plot_step'] == step]
        
        sns.boxplot(
            data=plot_df, 
            x='depth', 
            y='value', 
            hue='method', 
            order=depth_order, 
            palette=palette, 
            ax=ax,
            flierprops={"marker": "o", "markersize": 5, "alpha": 0.3}
            )
        
        y_label = f"Avg {resource_type.capitalize()} Usage"
        if resource_type == 'time':
            y_label += " (sec)"
            ax.set_yscale('log')
        elif resource_type == 'memory':
            y_label += " (kb)"
            ax.set_yscale('log')
        elif resource_type == 'cpu':
            y_label += " (%)"
        
        ax.set_title(f"{panel_labels[i]}) {step}", loc='left', fontsize=14)
        ax.set_ylabel(y_label)
        ax.set_xlabel("Depth Target" if i >= num_steps - num_cols else "")
        ax.grid(True, linestyle='--', alpha=0.6)
        ax.get_legend().remove()

    for i in range(num_steps, len(axes)):
        axes[i].set_visible(False)

    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, title='Method', loc='upper right', bbox_to_anchor=(0.99, 0.98))

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    filename = f"Supplemental_Figure_{resource_type.capitalize()}_Usage.pdf"
    plt.savefig(os.path.join(output_dir, filename), dpi=300, bbox_inches="tight", format="pdf")
    plt.show()

# --- 3. Main Execution ---
if __name__ == "__main__":
    main_df = load_and_reshape_data("data/total_resources.csv")
    
    output_directory = "supplemental_figures"
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    print(f"Saving supplemental figures to '{output_directory}/'...")

    for resource in ['cpu', 'memory', 'time']:
        print(f"\nGenerating multi-panel figure for {resource.upper()}...")
        create_multi_panel_supplemental_figure(main_df, resource, output_directory)

    print("\nScript finished.")