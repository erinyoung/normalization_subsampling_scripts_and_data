import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
import re

def parse_col_name(col_name):
    """
    Parses a column name into its constituent parts by first standardizing
    seqkit read counts to their corresponding depths.
    """
    
    num_depth_mapping = {
        990166: 10000, 495083: 5000, 99017: 1000, 49508: 500,
        29705: 400, 19803: 300, 14852: 200, 9902: 100,
        4951: 50, 2970: 30, 990: 10
    }
    
    # Identify resource from the end
    parts = col_name.split('_')
    resource = parts[-1]
    if resource not in ['cpu', 'time', 'memory']:
        return None, None, None

    # Identify method
    method = None
    if 'bbnorm' in col_name:
        method = 'bbnorm'
    elif 'seqkit' in col_name:
        method = 'seqkit'
    elif 'raw' in col_name:
        method = 'raw'
    else:
        return None, None, None

    # Handle the simple 'raw' case
    if method == 'raw':
        return 'raw', 'all', resource

    # --- Standardization Step ---
    # Create a temporary version of the column name to work with.
    # Replace any read count numbers with their corresponding depth numbers.
    standardized_col = col_name
    if method == 'seqkit':
        for read_count, depth_val in num_depth_mapping.items():
            if str(read_count) in standardized_col:
                standardized_col = standardized_col.replace(str(read_count), str(depth_val))

    # --- Extraction Step ---
    # Now that the column is standardized, extract the single numeric depth value.
    depth = None
    for part in standardized_col.split('_'):
        if part.isdigit():
            depth = part
            break
            
    if depth is None:
        return None, None, None # Unparsable column
            
    return method, str(depth), resource

def parse_time_to_seconds(time_str):
    """Converts a 'M:SS.ss' or 'S.ss' string to total seconds."""
    if pd.isna(time_str):
        return None
    if isinstance(time_str, (int, float)):
        return float(time_str)
        
    parts = str(time_str).split(':')
    if len(parts) == 2:
        minutes, seconds = parts
        return int(minutes) * 60 + float(seconds)
    elif len(parts) == 1:
        return float(parts[0])
    return None

def preprocess_resource_data(filepath):
    """
    Loads resource data, aggregates it by pipeline run (method + depth), and 
    transforms it into a long format suitable for analysis. Calculates sum for time,
    and both sum and mean for CPU and memory.

    Args:
        filepath (str): The path to the input CSV file ('total_resources.csv').

    Returns:
        pandas.DataFrame: A DataFrame where each row represents the total and mean 
                          resource usage for a single sample and pipeline combination.
    """
    df_wide = pd.read_csv(filepath)
    
    # 1. Melt the DataFrame from wide to long format
    df_long = pd.melt(df_wide, id_vars=['sample'], var_name='variable', value_name='value')
    df_long.dropna(subset=['value'], inplace=True)

    # 2. Parse the 'variable' column to extract method, depth, and metric
    df_long[['method', 'depth', 'metric']] = df_long['variable'].apply(
        lambda x: pd.Series(parse_col_name(x))
    )
    
    # Drop any rows where parsing failed
    df_long.dropna(subset=['method', 'depth', 'metric'], inplace=True)
    df_long.drop(columns=['variable'], inplace=True)

    # 3. Convert values to consistent numeric types
    time_mask = df_long['metric'] == 'time'
    df_long.loc[time_mask, 'value'] = df_long.loc[time_mask, 'value'].apply(parse_time_to_seconds)
    
    df_long['value'] = pd.to_numeric(df_long['value'], errors='coerce')
    df_long.dropna(subset=['value'], inplace=True)

    # 4. Aggregate (sum/mean) the values for each pipeline run
    # Separate data by metric to apply different aggregations
    df_cpu = df_long[df_long['metric'] == 'cpu']
    df_memory = df_long[df_long['metric'] == 'memory']
    df_time = df_long[df_long['metric'] == 'time']

    # Aggregate each metric type
    agg_cpu = df_cpu.groupby(['sample', 'method', 'depth'])['value'].agg(['sum', 'mean']).reset_index()
    agg_memory = df_memory.groupby(['sample', 'method', 'depth'])['value'].agg(['sum', 'mean']).reset_index()
    agg_time = df_time.groupby(['sample', 'method', 'depth'])['value'].agg('sum').reset_index()

    # Rename columns for clarity before merging
    agg_cpu.rename(columns={'sum': 'sum_cpu', 'mean': 'mean_cpu'}, inplace=True)
    agg_memory.rename(columns={'sum': 'sum_memory', 'mean': 'mean_memory'}, inplace=True)
    agg_time.rename(columns={'value': 'sum_time_seconds'}, inplace=True)

    # 5. Merge the aggregated dataframes back together
    df_final = pd.merge(agg_cpu, agg_memory, on=['sample', 'method', 'depth'], how='outer')
    df_final = pd.merge(df_final, agg_time, on=['sample', 'method', 'depth'], how='outer')
    
    return df_final


def compare_resources_to_control(df):
    """
    Compares resource usage of each method/depth group to a control ('raw'/'all')
    using t-tests and applies FDR correction.

    Args:
        df (pandas.DataFrame): Preprocessed DataFrame from preprocess_resource_data.

    Returns:
        pandas.DataFrame: A DataFrame containing the statistical comparison results,
                          including raw and adjusted p-values for all metrics.
    """
    # Define the control group
    control_group = df[(df['method'] == 'raw') & (df['depth'] == 'all')]
    
    # Get the data for the control group, dropping NaNs for comparison
    control_sum_cpu = control_group['sum_cpu'].dropna().astype(float)
    control_mean_cpu = control_group['mean_cpu'].dropna().astype(float)
    control_sum_time = control_group['sum_time_seconds'].dropna().astype(float)
    control_sum_memory = control_group['sum_memory'].dropna().astype(float)
    control_mean_memory = control_group['mean_memory'].dropna().astype(float)

    # Identify all other groups to test against the control
    test_groups = df[~((df['method'] == 'raw') & (df['depth'] == 'all'))]
    unique_groups = test_groups[['method', 'depth']].drop_duplicates()

    results = []
    
    # Perform t-test for each group against the control
    for index, row in unique_groups.iterrows():
        method, depth = row['method'], row['depth']
        current_group = test_groups[(test_groups['method'] == method) & (test_groups['depth'] == depth)]
        
        # Get data for the current group
        current_sum_cpu = current_group['sum_cpu'].dropna().astype(float)
        current_mean_cpu = current_group['mean_cpu'].dropna().astype(float)
        current_sum_time = current_group['sum_time_seconds'].dropna().astype(float)
        current_sum_memory = current_group['sum_memory'].dropna().astype(float)
        current_mean_memory = current_group['mean_memory'].dropna().astype(float)

        # Perform t-tests, handling cases with insufficient data
        p_sum_cpu = ttest_ind(control_sum_cpu, current_sum_cpu, nan_policy='omit', equal_var=False).pvalue if len(current_sum_cpu) > 1 and len(control_sum_cpu) > 1 else 1.0
        p_mean_cpu = ttest_ind(control_mean_cpu, current_mean_cpu, nan_policy='omit', equal_var=False).pvalue if len(current_mean_cpu) > 1 and len(control_mean_cpu) > 1 else 1.0
        p_sum_time = ttest_ind(control_sum_time, current_sum_time, nan_policy='omit', equal_var=False).pvalue if len(current_sum_time) > 1 and len(control_sum_time) > 1 else 1.0
        p_sum_memory = ttest_ind(control_sum_memory, current_sum_memory, nan_policy='omit', equal_var=False).pvalue if len(current_sum_memory) > 1 and len(control_sum_memory) > 1 else 1.0
        p_mean_memory = ttest_ind(control_mean_memory, current_mean_memory, nan_policy='omit', equal_var=False).pvalue if len(current_mean_memory) > 1 and len(control_mean_memory) > 1 else 1.0

        results.append({
            'method': method, 'depth': depth,
            'pvalue_sum_cpu': p_sum_cpu, 'pvalue_mean_cpu': p_mean_cpu,
            'pvalue_sum_time': p_sum_time,
            'pvalue_sum_memory': p_sum_memory, 'pvalue_mean_memory': p_mean_memory
        })

    if not results:
        return pd.DataFrame()

    results_df = pd.DataFrame(results)

    # --- Multiple Test Correction ---
    all_pvalues = (
        results_df['pvalue_sum_cpu'].tolist() + results_df['pvalue_mean_cpu'].tolist() +
        results_df['pvalue_sum_time'].tolist() +
        results_df['pvalue_sum_memory'].tolist() + results_df['pvalue_mean_memory'].tolist()
    )
    
    reject, pvals_corrected, _, _ = multipletests(all_pvalues, alpha=0.05, method='fdr_bh')
    
    num_comparisons = len(results_df)
    results_df['pvalue_adj_sum_cpu'] = pvals_corrected[0*num_comparisons : 1*num_comparisons]
    results_df['pvalue_adj_mean_cpu'] = pvals_corrected[1*num_comparisons : 2*num_comparisons]
    results_df['pvalue_adj_sum_time'] = pvals_corrected[2*num_comparisons : 3*num_comparisons]
    results_df['pvalue_adj_sum_memory'] = pvals_corrected[3*num_comparisons : 4*num_comparisons]
    results_df['pvalue_adj_mean_memory'] = pvals_corrected[4*num_comparisons : 5*num_comparisons]
        
    return results_df

if __name__ == '__main__':
    # Define input and output file paths
    input_filepath = 'data/total_resources.csv'
    output_filepath = 'resource_stats_results.csv'

    # Step 1: Load and preprocess the data
    print(f"Loading and preprocessing data from '{input_filepath}'...")
    preprocessed_df = preprocess_resource_data(input_filepath)
    print("Preprocessing complete.")
    
    # Step 2: Perform statistical analysis
    print("Performing statistical comparisons against the 'raw/all' control group...")
    stats_results_df = compare_resources_to_control(preprocessed_df)
    print("Analysis complete.")

    # Step 3: Save and display results
    stats_results_df.to_csv(output_filepath, index=False)
    print(f"\nResults saved to '{output_filepath}'")
    
    # Display results sorted by significance for time usage
    print("\n--- Statistical Analysis Results (sorted by adjusted p-value for total time) ---")
    print(stats_results_df.sort_values(by='pvalue_adj_sum_time').to_string())

