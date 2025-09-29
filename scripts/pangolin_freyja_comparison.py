import pandas as pd
import io

freyja_df = pd.read_csv("data/freyja_results.csv")

pango_df = pd.read_csv("data/pango_lineages.csv")

# Filter the pango_df for method=='all' and depth=='all'
pango_df_filtered = pango_df[(pango_df['method'] == 'raw') & (pango_df['depth'] == 'all')]
print(pango_df_filtered)

    # Merge the two dataframes
merged_df = pd.merge(freyja_df, pango_df_filtered, left_on='Sample', right_on='sample', how='inner')

# Compare the lineages
merged_df['lineage_match'] = merged_df['primary lineage'] == merged_df['lineage']

# Display the results
print("Comparison of Freyja and Pangolin Lineages (method=='all' and depth=='all'):")
print(merged_df[['Sample', 'primary lineage', 'lineage', 'method', 'depth', 'lineage_match']])

    # Display a summary
print("\nSummary:")
print(merged_df['lineage_match'].value_counts())