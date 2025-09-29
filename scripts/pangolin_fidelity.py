import pandas as pd

# This script requires the 'data/pango_lineages.csv' file to be present.
df = pd.read_csv("data/pango_lineages.csv")

# 1. Isolate the control data (method='raw' and depth='all')
control_df = df[(df['method'] == 'raw') & (df['depth'] == 'all')].copy()

# 2. Select and rename columns for clarity
control_df = control_df[['sample', 'lineage', 'pango_lineage_unaliased']].rename(
    columns={
        'lineage': 'control_lineage',
        'pango_lineage_unaliased': 'control_unaliased_lineage'
    }
)

# 3. Merge the control lineages back into the original dataframe
df_merged = pd.merge(df, control_df, on='sample', how='left')

# 4. Filter out the control rows to avoid self-comparison
comparison_df = df_merged[~((df_merged['method'] == 'raw') & (df_merged['depth'] == 'all'))].copy()


# 5. Define the new classification function with parental logic
def classify_mismatch(row):
    """
    Classifies the mismatch type with four categories:
    1. No Mismatch: Aliased lineages are identical.
    2. Aliased Mismatch: Aliased differ, but unaliased are identical.
    3. Parental Mismatch: Unaliased have a parent/sublineage relationship.
    4. True Mismatch: Unaliased are different and unrelated.
    """
    # Handle potential missing data by converting to string
    lineage = str(row['lineage'])
    control_lineage = str(row['control_lineage'])
    unaliased = str(row['pango_lineage_unaliased'])
    control_unaliased = str(row['control_unaliased_lineage'])

    # Category 1: Perfect match on the primary (aliased) lineage call
    if lineage == control_lineage:
        return 'No Mismatch'

    if lineage == "Unassigned":
        return 'Classification Failure'

    if control_lineage == "Unassigned":
        return 'unassigned to assigned'

    # Category 2: Aliased lineages differ, but unaliased are identical
    if unaliased == control_unaliased:
        return 'Aliased Mismatch'

    # Category 3: Check for parental relationship on unaliased lineages
    # We add a '.' to ensure we match whole numbers (e.g., B.1 is parent of B.1.1, not B.11)
    if control_unaliased.startswith(unaliased + '.'):
        return 'Gain of Specificity'

    if unaliased.startswith(control_unaliased + '.'):
        return 'Loss of Specificity'

    # Category 4: If none of the above, it's a significant, non-parental difference
    return 'Other Mismatch'

# 6. Apply the new classification logic
comparison_df['mismatch_type'] = comparison_df.apply(classify_mismatch, axis=1)

# 7. Display the rows that have any kind of mismatch
mismatched_rows = comparison_df[comparison_df['mismatch_type'] != 'No Mismatch']

print("--- Mismatched Lineage Report ---")
if not mismatched_rows.empty:
    print(mismatched_rows[[
        'sample', 'method', 'depth', 'lineage', 'control_lineage',
        'pango_lineage_unaliased', 'control_unaliased_lineage', 'mismatch_type'
    ]])
else:
    print("No mismatches found.")


# 8. Generate the detailed summary report with the new category
print("\n" + "="*65)
print("                  Detailed Mismatch Summary Report")
print("="*65)

summary_table = pd.crosstab(
    index=[comparison_df['method'], comparison_df['depth']],
    columns=comparison_df['mismatch_type']
)

if summary_table.empty:
    print("No comparisons were made.")
else:
    # Define the desired column order for the final report
    final_columns = ['No Mismatch', 'Loss of Specificity', 'Gain of Specificity', 'Other Mismatch', 'Classification Failure']
    
    # Ensure all category columns exist, even if there are no instances
    for col in final_columns:
        if col not in summary_table.columns:
            summary_table[col] = 0
    
    # Reorder columns for logical presentation
    print(summary_table[final_columns])



import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# This code assumes you have a pandas DataFrame named 'summary_table'
# with a multi-index of ('method', 'depth').

# --- 1. Prepare the DataFrame for Plotting ---
# Reset the index to make 'method' and 'depth' regular columns
df_for_plotting = summary_table.reset_index()
# Convert the 'depth' column to a numeric type
df_for_plotting['depth'] = pd.to_numeric(df_for_plotting['depth'])
# Sort the values correctly
df_for_plotting = df_for_plotting.sort_values(by=['method', 'depth'])
# Set the index back for plotting
df_sorted = df_for_plotting.set_index(['method', 'depth'])


# Define a logical order for the stacks (from bottom to top).
plot_columns = [
    'No Mismatch',
    'Loss of Specificity',
    'Gain of Specificity',
    'Other Mismatch',
    'Classification Failure'
]
# Ensure all columns exist in the DataFrame.
for col in plot_columns:
    if col not in df_sorted.columns:
        df_sorted[col] = 0

# --- 2. Create the Plot ---
fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(12, 10), sharey=True)
fig.suptitle('Lineage Classification Outcomes by Method and Depth', fontsize=16, y=0.98)

methods = df_sorted.index.get_level_values('method').unique()
colors = plt.cm.viridis(np.linspace(0, 1, len(plot_columns)))

for i, method in enumerate(methods):
    ax = axes[i]
    # Filter the data for the current method using .loc
    method_df = df_sorted.loc[method]
    
    method_df[plot_columns].plot(
        kind='bar',
        stacked=True,
        ax=ax,
        color=colors,
        width=0.8
    )
    
    ax.set_title(f'Method: {method}', fontsize=12)
    ax.set_ylabel('Number of Samples')
    ax.set_xlabel('Sequencing Depth')
    ax.tick_params(axis='x', rotation=45)
    ax.grid(axis='y', linestyle='--', alpha=0.7)
    
    if ax.get_legend():
        ax.get_legend().remove()

# --- 3. Add a Shared Legend and Save the Figure ---
handles, labels = ax.get_legend_handles_labels()
fig.legend(handles, labels, title='Mismatch Type', bbox_to_anchor=(1.05, 0.7))

plt.tight_layout(rect=[0, 0, 0.9, 0.96])

plt.savefig(
    'classification_summary.png',
    dpi=300,
    bbox_inches='tight'
)

print("Chart has been saved as 'classification_summary.png'")


import pandas as pd
from scipy.stats import chi2_contingency
from statsmodels.stats.multitest import multipletests
import numpy as np

# This code assumes you have your data in a DataFrame named 'df_sorted'
# with a multi-index of ('method', 'depth').

# --- Statistical Analysis ---

# 1. Set up the "perfect" control group robustly.
total_samples = df_sorted.iloc[0].sum()

# --- FIX STARTS HERE ---
# Create a Series of all zeros with the correct index
control_row = pd.Series(0, index=df_sorted.columns)
# Explicitly set the count in the 'No Mismatch' column, regardless of its position
control_row['No Mismatch'] = total_samples
# --- FIX ENDS HERE ---


# 2. Loop through each experimental group and perform the test.
results = []
for index, row in df_sorted.iterrows():
    # Combine rows and filter out uninformative zero-sum columns
    temp_table = pd.DataFrame([row, control_row])
    zero_sum_cols = temp_table.columns[temp_table.sum() == 0]
    contingency_table = temp_table.drop(columns=zero_sum_cols).values

    # Ensure the table is still valid before running the test
    if contingency_table.shape[1] > 1:
        chi2, p_value, _, _ = chi2_contingency(contingency_table)
    else:
        p_value = 1.0

    results.append({
        'method': index[0],
        'depth': index[1],
        'p_value': p_value
    })

results_df = pd.DataFrame(results)

# 3. Apply the Bonferroni correction for multiple comparisons.
reject, p_values_corrected, _, _ = multipletests(
    results_df['p_value'],
    alpha=0.05,
    method='bonferroni'
)

results_df['p_value_corrected'] = p_values_corrected
results_df['is_significant'] = results_df['p_value_corrected'] < 0.05

# Sort for clarity
results_df = results_df.sort_values(by=['method', 'depth'])

print("Statistical Significance Results:")
print(results_df)
