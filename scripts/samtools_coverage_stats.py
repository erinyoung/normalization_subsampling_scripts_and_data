import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

def compare_to_control(df):
    """
    Compare method/depth combinations to a control group using t-tests and FDR correction.

    This function compares coverage-related metrics between each method/depth 
    combination in the input DataFrame and a designated control group 
    (method = 'all', depth = 'all'). Independent two-sample t-tests are used, 
    and p-values are adjusted for multiple testing using the 
    Benjamini-Hochberg false discovery rate (FDR) procedure.

    Parameters
    ----------
    df : pandas.DataFrame
        Input DataFrame containing the following columns:
        - 'method' (str): Method label.
        - 'depth' (str): Depth label.
        - 'covbases' (numeric): Number of covered bases.
        - 'coverage' (numeric): Coverage fraction or percentage.
        - 'meandepth' (numeric): Mean sequencing depth.
        - 'meanbaseq' (numeric): Mean base quality.
        - 'meanmapq' (numeric): Mean mapping quality.

    Returns
    -------
    pandas.DataFrame
        A DataFrame with one row per method/depth combination (excluding the control).
        Includes raw and FDR-adjusted p-values for each metric:
        - pvalue_covbases, pvalue_coverage, pvalue_meandepth, pvalue_meanbaseq, pvalue_meanmapq
        - covbases_pvalue_adj, coverage_pvalue_adj, meandepth_pvalue_adj, meanbaseq_pvalue_adj, meanmapq_pvalue_adj

    Notes
    -----
    - Control group is defined as rows where method == 'all' and depth == 'all'.
    - Independent two-sample t-tests are performed with `nan_policy='omit'`.
    - Multiple test correction is applied across all metrics and comparisons.
    """

    #covbases,coverage,meandepth,meanbaseq,meanmapq
    
    # Define the control groups
    control_group     = df[(df['method'] == 'all') & (df['depth'] == 'all')]
    control_covbases  = control_group['covbases']
    control_coverage  = control_group['coverage']
    control_meandepth = control_group['meandepth']
    control_meanbaseq = control_group['meanbaseq']
    control_meanmapq  = control_group['meanmapq']
    
    # Get unique method and depth combinations, excluding the control
    combinations = df[['method', 'depth']].drop_duplicates()
    combinations = combinations[~((combinations['method'] == 'all') & (combinations['depth'] == 'all'))]

    results = []

    # Iterate over each combination and perform t-test
    for index, row in combinations.iterrows():
        method = row['method']
        depth  = row['depth']
        
        current_group = df[(df['method'] == method) & (df['depth'] == depth)]
        
        ttest_covbases  = ttest_ind(control_covbases, current_group['covbases'], nan_policy='omit')    
        ttest_coverage  = ttest_ind(control_coverage, current_group['coverage'], nan_policy='omit')
        ttest_meandepth = ttest_ind(control_meandepth, current_group['meandepth'], nan_policy='omit')
        ttest_meanbaseq = ttest_ind(control_meanbaseq, current_group['meanbaseq'], nan_policy='omit')
        ttest_meanmapq  = ttest_ind(control_meanmapq, current_group['meanmapq'], nan_policy='omit')

        results.append({
            'method'           : method,
            'depth'            : depth,
            'pvalue_covbases'  : ttest_covbases.pvalue,
            'pvalue_coverage'  : ttest_coverage.pvalue,
            'pvalue_meandepth' : ttest_meandepth.pvalue,
            'pvalue_meanbaseq' : ttest_meanbaseq.pvalue,
            'pvalue_meanmapq'  : ttest_meanmapq.pvalue
        })

    results_df = pd.DataFrame(results)

    # --- Multiple Test Correction ---
    # Combine all p-values into one list for correction
    all_pvalues = results_df['pvalue_covbases'].tolist() + results_df['pvalue_coverage'].tolist() + results_df['pvalue_meandepth'].tolist() + results_df['pvalue_meanbaseq'].tolist() + results_df['pvalue_meanmapq'].tolist()
    
    # Apply Benjamini-Hochberg FDR correction
    reject, pvals_corrected, _, _ = multipletests(all_pvalues, alpha=0.05, method='fdr_bh')
    
    # Split the corrected p-values back into their respective columns
    num_comparisons = len(results_df)
    results_df['covbases_pvalue_adj']  = pvals_corrected[0*num_comparisons : 1*num_comparisons]
    results_df['coverage_pvalue_adj']  = pvals_corrected[1*num_comparisons : 2*num_comparisons]
    results_df['meandepth_pvalue_adj'] = pvals_corrected[2*num_comparisons : 3*num_comparisons]
    results_df['meanbaseq_pvalue_adj'] = pvals_corrected[3*num_comparisons : 4*num_comparisons]
    results_df['meanmapq_pvalue_adj']  = pvals_corrected[4*num_comparisons : 5*num_comparisons]
        
    return results_df

# Load the dataset
df = pd.read_csv("data/samtools_coverage.csv")

pd.options.display.float_format = '{:.4f}'.format

# Get the results
results_df = compare_to_control(df)

# Print the results
print("Comparison of method/depth combinations to control with FDR correction:")
print(results_df)

# Save the results to a CSV file
results_df.to_csv('samtools_coverage_pvalues_corrected.csv', index=False)