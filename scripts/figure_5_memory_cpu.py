import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

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

# --- 2. Plotting Function for Figure 5 ---

def plot_resource_boxplots(df):
    """Creates the side-by-side boxplots for Memory and CPU (Figure 5)."""
    fig, axes = plt.subplots(1, 2, figsize=(22, 8))
    fig.suptitle("Computational Resource Utilization by Method and Depth", fontsize=18, y=0.98)

    palette = {'raw': '#0072B2', 'bbnorm': '#D55E00', 'seqkit': '#009E73'}
    depth_order = ['10', '30', '50', '100', '200', '300', '400', '500', '1000', '5000', '10000', 'all']

    # --- Panel A: Memory Usage ---
    df_mem = df[(df['resource'] == 'memory') & (df['depth'].isin(depth_order))].copy()
    df_mem['method'] = df_mem['method'].replace({'bbnorm': 'bbnorm', 'seqkit': 'seqkit', 'raw': 'raw'})
    
    sns.boxplot(
        data=df_mem, 
        x='depth', 
        y='value', 
        hue='method', 
        order=depth_order, 
        palette=palette, 
        ax=axes[0],
        flierprops={"marker": "o", "markersize": 5, "alpha": 0.3}
        )
    axes[0].set_title("A) Distribution of Memory Usage", loc='left', fontsize=16)
    axes[0].set_ylabel("Average Memory Usage (kbytes)")
    axes[0].set_xlabel("Subsampling Depth Target")
    axes[0].set_yscale('log')
    axes[0].grid(True, linestyle='--', alpha=0.6)

    # --- Panel B: CPU Usage ---
    df_cpu = df[(df['resource'] == 'cpu') & (df['depth'].isin(depth_order))].copy()
    df_cpu['method'] = df_cpu['method'].replace({'bbnorm': 'bbnorm', 'seqkit': 'seqkit', 'raw': 'raw'})

    sns.boxplot(
        data=df_cpu, 
        x='depth',
        y='value', 
        hue='method', 
        order=depth_order, 
        palette=palette, 
        ax=axes[1],
        flierprops={"marker": "o", "markersize": 5, "alpha": 0.3}
        )
    axes[1].set_title("B) Distribution of CPU Usage", loc='left', fontsize=16)
    axes[1].set_ylabel("Average CPU Usage (percent)")
    axes[1].set_xlabel("Subsampling Depth Target")
    axes[1].grid(True, linestyle='--', alpha=0.6)

    # --- Final Touches: Shared Legend ---
    # Get handles from one plot and remove individual legends
    handles, labels = axes[0].get_legend_handles_labels()
    axes[0].get_legend().remove()
    axes[1].get_legend().remove()
    
    # Add a single, shared legend, anchored to the right of Panel B
    axes[1].legend(handles, labels, title='Method', loc='upper left', bbox_to_anchor=(1.02, 1))
    
    # Adjust layout to make space for the legend
    plt.tight_layout(rect=[0, 0, 0.88, 0.93])
    plt.savefig("Figure5.tiff", dpi=300, bbox_inches="tight", format="tiff", pil_kwargs={"compression": "tiff_lzw"})

# --- 3. Main Execution ---
if __name__ == "__main__":
    main_df = load_and_reshape_data("data/total_resources.csv")
    
    print("Generating Figure 5: Resource Utilization Boxplots...")
    plot_resource_boxplots(main_df)
    
    print("\nScript finished.")
