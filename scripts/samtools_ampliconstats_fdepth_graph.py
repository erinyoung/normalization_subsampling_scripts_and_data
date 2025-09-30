import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_csv("data/ampliconstats_fdepth.csv")

# Identify the amplicon columns
amplicon_columns = df.columns[df.columns.get_loc('Depth') + 1:]

# Calculate counts for each sample (row)
df['Failed'] = (df[amplicon_columns] < 10).sum(axis=1)
df['Depth >=10 and <=50'] = ((df[amplicon_columns] >= 10) & (df[amplicon_columns] <= 50)).sum(axis=1)
df['Passed'] = (df[amplicon_columns] > 50).sum(axis=1)

# Group by Method and Depth
grouped = df.groupby(['Method', 'Depth'])

# Calculate the mean and standard deviation for our categories
agg_df = grouped[['Failed', 'Depth >=10 and <=50', 'Passed']].agg(['mean', 'std'])

# Flatten the multi-index columns (e.g., ('Passed', 'mean') -> 'Passed_mean')
agg_df.columns = ['_'.join(col).strip() for col in agg_df.columns.values]

# Fill NaN in stddev columns with 0 (for groups with only one sample)
agg_df = agg_df.fillna(0).reset_index()
#agg_df = agg_df.sort_values(by='Depth')
agg_df['Depth'] = agg_df['Depth'].replace('all', np.inf)
agg_df['Depth'] = pd.to_numeric(agg_df['Depth'])
agg_df = agg_df.sort_values(by=['Method', 'Depth'])
print(agg_df)

agg_df['label'] = agg_df['Method'] + "_" + agg_df['Depth'].astype(str).str.replace("raw_inf","raw").str.replace('.0', '', regex=False)

agg_df = agg_df.drop(columns=['Depth', 'Method', 'Depth >=10 and <=50_std', 'Failed_std', 'Passed_std'])

agg_df.set_index('label', inplace=True)

agg_df.rename(columns={
    'Failed_mean': 'Failed',
    'Depth >=10 and <=50_mean': 'Depth >=10 and <=50',
    'Passed_mean': 'Passed'
}, inplace=True)


agg_df = agg_df[agg_df.columns[::-1]]

ax = agg_df.plot(kind='bar',
            stacked=True,
            figsize=(12, 8),
            edgecolor='black',
            linewidth=1,
            color=['black', 'gray', 'white'])

# 4. Customize the plot
ax.set_title('Amplicon Pass/Fail for each Method and Depth', fontsize=16)
ax.set_ylabel('Average Number of Amplicons')
ax.set_xlabel('Method and Depth')
ax.legend(title='Depth Category', bbox_to_anchor=(1.02, 1), loc='upper left')

# Ensure the plot layout is clean
#plt.tight_layout()
plt.savefig(
    'samtools_ampliconstats_fdepth_counts_average_plot.png', 
    bbox_inches='tight', 
    dpi=300
    )
plt.show()
