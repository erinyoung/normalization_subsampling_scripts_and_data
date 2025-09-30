import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

def compare_to_control(df):
    """
    Compare scorpio metrics of method/depth combinations to a control group with FDR correction.

    This function compares three metrics (`scorpio_support`, `scorpio_conflict`, 
    and extracted `ambiguity_content`) between each method/depth combination and 
    a designated control group (method = 'raw', depth = 'all'). Independent 
    two-sample t-tests are performed, and p-values are adjusted for multiple 
    testing using the Benjamini-Hochberg false discovery rate (FDR) procedure.

    Parameters
    ----------
    df : pandas.DataFrame
        Input DataFrame containing at least the following columns:
        - 'method' (str): Method label.
        - 'depth' (str): Depth label.
        - 'scorpio_support' (numeric): Scorpio support values.
        - 'scorpio_conflict' (numeric): Scorpio conflict values.
        - 'qc_notes' (str): QC notes string, from which ambiguity values (formatted as 'Ambiguous_content:<float>') will be extracted.

    Returns
    -------
    pandas.DataFrame
        DataFrame with one row per method/depth combination (excluding the control).
        Includes raw and FDR-adjusted p-values for each metric:
        - scorpio_support_pvalue, scorpio_conflict_pvalue, ambiguity_content_pvalue
        - scorpio_support_pvalue_adj, scorpio_conflict_pvalue_adj, ambiguity_content_pvalue_adj

    Notes
    -----
    - Control group is defined as rows where method == 'raw' and depth == 'all'.
    - Ambiguity values are parsed from the `qc_notes` column using regex.
    - Independent two-sample t-tests are performed with `nan_policy='omit'`.
    - Multiple test correction is applied jointly across all metrics and comparisons.
    """


    df["ambiguity_content"] = df["qc_notes"].str.extract(r"Ambiguous_content:([0-9.]+)")[0].astype(float)

    
    # Define the control group
    control_group = df[(df['method'] == 'raw') & (df['depth'] == 'all')]
    control_support = control_group['scorpio_support']
    control_conflict = control_group['scorpio_conflict']
    control_ambiguity = control_group['ambiguity_content']
    
    # Get unique method and depth combinations, excluding the control
    combinations = df[['method', 'depth']].drop_duplicates()
    combinations = combinations[~((combinations['method'] == 'raw') & (combinations['depth'] == 'all'))]

    results = []

    # Iterate over each combination and perform t-test
    for index, row in combinations.iterrows():
        method = row['method']
        depth = row['depth']
        
        current_group = df[(df['method'] == method) & (df['depth'] == depth)]
        
        # T-test for scorpio_support
        support_ttest = ttest_ind(control_support, current_group['scorpio_support'], nan_policy='omit')
        
        # T-test for scorpio_conflict
        conflict_ttest = ttest_ind(control_conflict, current_group['scorpio_conflict'], nan_policy='omit')

        # T-test for ambiguity
        ambiguity_ttest = ttest_ind(control_ambiguity, current_group['ambiguity_content'], nan_policy='omit')
        
        results.append({
            'method': method,
            'depth': depth,
            'scorpio_support_pvalue': support_ttest.pvalue,
            'scorpio_conflict_pvalue': conflict_ttest.pvalue,
            'ambiguity_content_pvalue': ambiguity_ttest.pvalue,
        })

    results_df = pd.DataFrame(results)

    # --- Multiple Test Correction ---
    # Combine all p-values into one list for correction
    all_pvalues = results_df['scorpio_support_pvalue'].tolist() + results_df['scorpio_conflict_pvalue'].tolist() +results_df['ambiguity_content_pvalue'].tolist() 
    
    # Apply Benjamini-Hochberg FDR correction
    reject, pvals_corrected, _, _ = multipletests(all_pvalues, alpha=0.05, method='fdr_bh')
    
    # Split the corrected p-values back into their respective columns
    num_comparisons = len(results_df)
    results_df['scorpio_support_pvalue_adj'] = pvals_corrected[:num_comparisons]
    results_df['scorpio_conflict_pvalue_adj'] = pvals_corrected[num_comparisons:2*num_comparisons]
    results_df['ambiguity_content_pvalue_adj'] = pvals_corrected[2*num_comparisons:]
        
    return results_df

# Load the dataset
# Make sure 'pango_lineages.csv' is in the same directory
df = pd.read_csv('data/pango_lineages.csv')

pd.options.display.float_format = '{:.4f}'.format

# Get the results
results_df = compare_to_control(df)

# Print the results
print("Comparison of method/depth combinations to control with FDR correction:")
print(results_df)

# Save the results to a CSV file
results_df.to_csv('pangolin_comparison_pvalues_corrected.csv', index=False)