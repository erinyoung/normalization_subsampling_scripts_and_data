import pandas as pd
import numpy as np

# 1. Load the data
file_path = 'data/wastewater_freyja.tsv'
df = pd.read_csv(file_path, sep='\t')

# 2. Helper function to parse filenames
def parse_sample_name(raw_name):
    raw_name = raw_name.replace('_variants.tsv', '')
    if '_bbnorm_' in raw_name:
        base, depth = raw_name.split('_bbnorm_')
        return base, 'BBNorm', int(depth)
    elif '_seqkit_' in raw_name:
        base, depth = raw_name.split('_seqkit_')
        return base, 'SeqKit', int(depth)
    else:
        return raw_name, 'Original', float('inf')

df['BaseSample'], df['Method'], df['Depth'] = zip(*df['Unnamed: 0'].apply(parse_sample_name))

# 3. Establish the "ground truth" dominant lineage for each sample from the Original run
orig_dom = {}
for sample in df['BaseSample'].unique():
    orig_row = df[(df['BaseSample'] == sample) & (df['Method'] == 'Original')]
    if len(orig_row) > 0 and not pd.isna(orig_row.iloc[0]['lineages']):
        lins = str(orig_row.iloc[0]['lineages']).split()
        if lins:
            orig_dom[sample] = lins[0]

# 4. Aggregate metrics by Depth and Method
results = []
depths = [float('inf'), 10000, 5000, 1000, 500, 400, 300, 200, 100, 50, 30, 10]

for d in depths:
    methods_to_run = ['Original'] if d == float('inf') else ['BBNorm', 'SeqKit']
    
    for m in methods_to_run:
        subset = df[(df['Depth'] == d) & (df['Method'] == m)]
        depth_label = 'Original' if d == float('inf') else f"{int(d)}X"
        method_label = 'None' if m == 'Original' else m
        
        n_lins = []
        dom_abunds = []
        matches = 0
        total = 0
        
        for _, row in subset.iterrows():
            if pd.isna(row['lineages']): continue
            
            lins = str(row['lineages']).split()
            abunds = [float(x) for x in str(row['abundances']).split()]
            n_lins.append(len(lins))
            samp = row['BaseSample']
            
            # Track the original dominant lineage's survival
            if samp in orig_dom and orig_dom[samp] in lins:
                dom_abunds.append(abunds[lins.index(orig_dom[samp])])
                # Check if it is STILL the primary dominant lineage
                if lins[0] == orig_dom[samp]:
                    matches += 1
            else:
                dom_abunds.append(0.0)
            total += 1
            
        results.append({
            'Target Depth': depth_label,
            'Method': method_label,
            'Avg. Lineages Detected': round(np.mean(n_lins), 1) if n_lins else 0,
            'Original Dominant Abundance': f"{np.mean(dom_abunds)*100:.1f}%" if dom_abunds else "0.0%",
            'Dominant Match Rate': f"{(matches/total)*100:.0f}%" if total > 0 else "0%"
        })

# 5. Output the final DataFrame for Table 4
table_4_df = pd.DataFrame(results)
print(table_4_df.to_string(index=False))
table_4_df.to_csv('Table_4_Lineage_Shifts.csv', index=False)