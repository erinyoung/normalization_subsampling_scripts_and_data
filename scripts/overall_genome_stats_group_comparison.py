import pandas as pd
from scipy import stats
import scikit_posthocs as sp

# Load your dataframe
df = pd.read_csv("data/samtools_coverage.csv")

# Create a single column to identify each group
# Note: .astype(str) is used in case the depth column is numeric
df['group'] = df['method'] + '_' + df['depth'].astype(str)

for val_col in ["numreads","covbases","coverage","meandepth","meanbaseq","meanmapq"]:
    print(f"Getting p-values for {val_col}")

    # Perform Dunnett's test
    p_values = sp.posthoc_dunnett(df, val_col=val_col, group_col='group', control='all_all')

    print("Dunnett's Test P-Values (compared to control 'all_all'):")
    print(p_values["all_all"])
    p_values.to_csv(f"samtools_coverage_stats_{val_col}.csv", index=False)