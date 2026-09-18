import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

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

# --- 2. Plotting Function for Three-Panel Figure 4 ---

def plot_time_summary_figure_by_method(df):
    """Creates a three-panel figure for processing time, with one panel per method."""
    df_time = df[df['resource'] == 'time'].copy()
    
    step_order = [
        "bbnorm", "seqkit_R1", "seqkit_R2", "fastp", "bbmap", "sort",
        "ivar_trim", "sort2", "ivar_consensus_mpileup", "ivar_consensus_consensus"
    ]
    
    colors = cm.get_cmap('viridis', len(step_order))
    color_map = {step: colors(i) for i, step in enumerate(step_order)}
    
    avg_times = df_time.groupby(['method', 'depth', 'step'])['value'].mean().unstack().fillna(0)
    avg_times = avg_times.reset_index()
    
    depth_order_map = {d: i for i, d in enumerate(['10', '30', '50', '100', '200', '300', '400', '500', '1000', '5000', '10000', 'all'])}
    avg_times['depth_order'] = avg_times['depth'].map(depth_order_map)
    avg_times = avg_times.sort_values(by=['method', 'depth_order'])

    fig, axes = plt.subplots(3, 1, figsize=(18, 20), sharex=False)
    fig.suptitle("Average Pipeline Processing Time", fontsize=20, y=0.98)

    # --- Panel A: All Data ---
    all_df = avg_times.copy()
    all_df['method_depth'] = all_df['method'] + '_' + all_df['depth'].astype(str)
    all_df = all_df.set_index('method_depth')
    all_df[step_order].plot(kind='bar', stacked=True, ax=axes[0], width=0.8, legend=True, color=color_map)
    axes[0].set_title("A) Overall Comparison", loc='left', fontsize=16)
    axes[0].set_ylabel("Average Time (seconds)")
    # Position this panel's legend outside the plot
    axes[0].legend(title="Pipeline Step", bbox_to_anchor=(1.02, 1), loc='upper left')

    # --- Panel B: BBNorm Pipeline ---
    bbnorm_df = avg_times[avg_times['method'] == 'bbnorm'].set_index('depth').sort_values('depth_order')
    bbnorm_df[step_order].plot(kind='bar', stacked=True, ax=axes[1], width=0.8, legend=False, color=color_map)
    axes[1].set_title("B) BBNorm Pipeline", loc='left', fontsize=16)
    axes[1].set_ylabel("Average Time (seconds)")

    # --- Panel C: SeqKit Pipeline ---
    seqkit_df = avg_times[avg_times['method'] == 'seqkit'].set_index('depth').sort_values('depth_order')
    seqkit_df[step_order].plot(kind='bar', stacked=True, ax=axes[2], width=0.8, legend=False, color=color_map)
    axes[2].set_title("C) SeqKit Pipeline", loc='left', fontsize=16)
    axes[2].set_ylabel("Average Time (seconds)")
    axes[2].set_xlabel("Depth Target")
    
    # --- Final Touches ---
    plt.tight_layout(rect=[0, 0, 0.88, 0.95])
    plt.savefig("Figure4.tiff", dpi=300, format="tiff", pil_kwargs={"compression": "tiff_lzw"})

# --- 3. Main Execution ---
if __name__ == "__main__":
    main_df = load_and_reshape_data("data/total_resources.csv")
    
    print("Generating three-panel Figure 4 with corrected legend...")
    plot_time_summary_figure_by_method(main_df)
    
    print("\nScript finished.")