import pandas as pd
from scipy.stats import mannwhitneyu
import numpy as np

def analyze_coverage_evenness():
    """
    Analyzes the evenness of amplicon coverage for each Sample.

    This function calculates the standard deviation of amplicon proportions across
    all amplicons for each Sample as a measure of evenness. It then uses a
    Mann-Whitney U test to compare the evenness of each Method/Depth group
    against a baseline of all 'raw' data.

    Args:
        file_path (str): The path to the input CSV file.

    Returns:
        pandas.DataFrame: A DataFrame containing the results of the statistical tests, including the Method, Depth, test statistic, and p-value.
    """
    df = pd.read_csv("data/ampliconstats_frperc.csv")

    # Identify amplicon columns
    amplicon_cols = [col for col in df.columns if col.startswith("SARS-CoV-2")]
    
    print(f"Found {len(amplicon_cols)} amplicon columns to calculate evenness from.")

    # --- 2. Calculate Evenness Metric for Each Sample ---
    # We calculate the standard deviation *across the row* for the amplicon columns.
    # A smaller SD means more even coverage.
    df['evenness_sd'] = df[amplicon_cols].std(axis=1)
    
    print("\nCalculated 'evenness_sd' for each Sample.")
    print(df[['Sample', 'Method', 'Depth', 'evenness_sd']].head())

    # --- 3. Create Baseline and Compare Groups ---
    # Baseline group: 'raw' Method data
    baseline_evenness = df[df['Method'] == 'raw']['evenness_sd'].dropna()
    
    # Get unique Method/Depth combinations to test against the baseline
    test_groups = df[df['Method'] != 'raw'][['Method', 'Depth']].drop_duplicates()

    results = []
    for index, row in test_groups.iterrows():
        current_Method = row['Method']
        current_Depth = row['Depth']
        
        # Comparison group: data for the current Method and Depth
        comparison_evenness = df[
            (df['Method'] == current_Method) & (df['Depth'] == current_Depth)
        ]['evenness_sd'].dropna()
        
        # Ensure there's enough data to perform the test
        if len(comparison_evenness) < 1:
            continue
            
        # Perform Mann-Whitney U test. 'alternative="two-sided"' checks if the
        # distributions are different in general.
        stat, p_value = mannwhitneyu(
            baseline_evenness, 
            comparison_evenness, 
            alternative='two-sided'
        )
        
        results.append({
            'Method': current_Method,
            'Depth': current_Depth,
            'median_evenness_sd': comparison_evenness.median(),
            'baseline_median_sd': baseline_evenness.median(),
            'mwu_statistic': stat,
            'p_value': p_value
        })

    # --- 4. Format and Return Results ---
    results_df = pd.DataFrame(results)
    
    # Add a column to easily identify significant results (alpha = 0.05)
    results_df['is_significant'] = results_df['p_value'] < 0.05
    
    return results_df.sort_values('p_value')


evenness_results = analyze_coverage_evenness()

evenness_results.to_csv("samtools_ampliconstats_frperc_stats.csv")

if evenness_results is not None:
    print("\n--- Coverage Evenness Analysis Complete ---")
    print("A lower 'median_evenness_sd' indicates more even coverage.")
    print(evenness_results)
