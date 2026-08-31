import pandas as pd
import numpy as np
import re

# 1. Load data
file_path = 'data/wastewater_freyja.tsv'
df = pd.read_csv(file_path, sep='\t')

# 2. Helper functions
def get_base_lineage(lin):
    return str(lin).split('.')[0]

def parse_sample_info(raw_name):
    name = str(raw_name).replace('_variants.tsv', '').replace('.tsv', '')
    
    # Extract depth
    match = re.search(r'_(\d+)$', name)
    if match:
        depth = int(match.group(1))
        name = name[:match.start()]
    else:
        depth = 'Raw'
        
    # Extract method
    method = 'Raw'
    if 'bbnorm' in name.lower():
        method = 'BBNorm'
        name = re.sub(r'_?bbnorm', '', name, flags=re.IGNORECASE)
    elif 'seqkit' in name.lower():
        method = 'SeqKit'
        name = re.sub(r'_?seqkit', '', name, flags=re.IGNORECASE)
        
    base_sample = name.replace('-UT', '').rstrip('_')
    return base_sample, method, depth

# 3. First Pass: Identify Top 19 Overall Base Lineages
base_totals = {}
for idx, row in df.iterrows():
    if pd.isna(row['lineages']) or pd.isna(row['abundances']):
        continue
    lins = str(row['lineages']).split()
    abunds = str(row['abundances']).split()
    for l, a in zip(lins, abunds):
        try:
            v = float(a)
            b = get_base_lineage(l)
            base_totals[b] = base_totals.get(b, 0.0) + v
        except ValueError:
            pass

top_19_bases = set(sorted(base_totals, key=base_totals.get, reverse=True)[:19])

# 4. Parse all sample rows
parsed_records = []
for idx, row in df.iterrows():
    if pd.isna(row['lineages']) or pd.isna(row['abundances']):
        continue
        
    raw_name = row['Unnamed: 0']
    base_sample, method, depth = parse_sample_info(raw_name)
    
    lins = str(row['lineages']).split()
    abunds = str(row['abundances']).split()
    
    base_abundances = {}
    other_val = 0.0
    
    for l, a in zip(lins, abunds):
        try:
            v = float(a)
            b = get_base_lineage(l)
            base_abundances[b] = base_abundances.get(b, 0.0) + v
            if b not in top_19_bases:
                other_val += v
        except ValueError:
            pass
            
    dom_lineage = max(base_abundances, key=base_abundances.get) if base_abundances else 'None'
    
    parsed_records.append({
        'Base_Sample': base_sample,
        'Method': method,
        'Depth': depth,
        'Dominant_Lineage': dom_lineage,
        'Base_Abundances': base_abundances,
        'Other_Abundance': other_val
    })

rec_df = pd.DataFrame(parsed_records)

# 5. Extract Raw Baseline Reference for Comparison
raw_references = {}
raw_rows = rec_df[rec_df['Depth'] == 'Raw']
for idx, r in raw_rows.iterrows():
    raw_references[r['Base_Sample']] = {
        'dom': r['Dominant_Lineage'],
        'abundances': r['Base_Abundances']
    }

# 6. Compare Reduced Datasets against Raw References
comparison_results = []
for idx, r in rec_df.iterrows():
    # Skip raw rows themselves in comparison metrics
    if r['Depth'] == 'Raw':
        continue
        
    sample = r['Base_Sample']
    if sample not in raw_references:
        continue
        
    raw_ref = raw_references[sample]
    
    # Metric 1: Dominant Lineage Match (Concordance)
    is_match = (r['Dominant_Lineage'] == raw_ref['dom'])
    
    # Metric 2: Mean Lineage Abundance Difference (Delta)
    all_lineages = set(r['Base_Abundances'].keys()).union(set(raw_ref['abundances'].keys()))
    diffs = [abs(r['Base_Abundances'].get(k, 0.0) - raw_ref['abundances'].get(k, 0.0)) for k in all_lineages]
    mean_diff = np.mean(diffs) if diffs else 0.0
    
    # Metric 3: Other > 10%
    high_other = (r['Other_Abundance'] > 0.10)
    
    comparison_results.append({
        'Method': r['Method'],
        'Depth': r['Depth'],
        'Match': is_match,
        'Mean_Diff': mean_diff,
        'High_Other': high_other
    })

comp_df = pd.DataFrame(comparison_results)

# 7. Aggregate Summary Table
summary_table = comp_df.groupby(['Depth', 'Method']).agg(
    Concordance_Pct=('Match', lambda x: f"{np.mean(x) * 100:.1f}%"),
    Mean_Abund_Diff=('Mean_Diff', lambda x: f"{np.mean(x):.3f}"),
    High_Other_Count=('High_Other', 'sum')
).reset_index()

# Filter for key target depth benchmarks
target_depths = [5000, 1000, 500, 200, 100, 50, 10]
summary_table = summary_table[summary_table['Depth'].isin(target_depths)]

# Sort by Depth descending, then Method
summary_table['Depth'] = pd.Categorical(summary_table['Depth'], categories=target_depths, ordered=True)
summary_table = summary_table.sort_values(by=['Depth', 'Method'])

# Rename columns for manuscript publication
summary_table.columns = [
    'Target Depth', 
    'Reduction Method', 
    'Dominant Lineage Concordance (%)', 
    'Mean Lineage Abundance Difference (Δ)', 
    'Samples with "Other" > 10% (n)'
]

print("\n--- SUMMARY TABLE (OPTION 2) ---")
print(summary_table.to_string(index=False))

# Export to CSV for manuscript inclusion
summary_table.to_csv('table_2_summary_metrics.csv', index=False)