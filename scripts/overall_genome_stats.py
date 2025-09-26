import pandas as pd
import io
import scipy.stats as stats

df = pd.read_csv("data/samtools_coverage.csv")

#covbases,coverage,meandepth,meanbaseq,meanmapq

# Separate the covbases data for each group
df_all = df[df['method'] == 'all']['covbases']
df_bbnorm = df[df['method'] == 'bbnorm']['covbases']
df_seqkit = df[df['method'] == 'seqkit']['covbases']

print("="*40)
print("  Comparing covbases for Each Method to 'all'")
print("="*40)

# --- Test 1: 'all' vs 'bbnorm' ---
stat_bbnorm, p_bbnorm = stats.ttest_ind(df_all, df_bbnorm)
print(f"Comparison: 'all' vs 'bbnorm'")
print(f"P-value: {p_bbnorm:.4f}")

# --- Test 2: 'all' vs 'seqkit' ---
stat_seqkit, p_seqkit = stats.ttest_ind(df_all, df_seqkit)
print(f"\nComparison: 'all' vs 'seqkit'")
print(f"P-value: {p_seqkit:.4f}")

#covbases,coverage,meandepth,meanbaseq,meanmapq

# Separate the coverage data for each group
df_all = df[df['method'] == 'all']['coverage']
df_bbnorm = df[df['method'] == 'bbnorm']['coverage']
df_seqkit = df[df['method'] == 'seqkit']['coverage']

print("="*40)
print("  Comparing coverage for Each Method to 'all'")
print("="*40)

# --- Test 1: 'all' vs 'bbnorm' ---
stat_bbnorm, p_bbnorm = stats.ttest_ind(df_all, df_bbnorm)
print(f"Comparison: 'all' vs 'bbnorm'")
print(f"P-value: {p_bbnorm:.4f}")

# --- Test 2: 'all' vs 'seqkit' ---
stat_seqkit, p_seqkit = stats.ttest_ind(df_all, df_seqkit)
print(f"\nComparison: 'all' vs 'seqkit'")
print(f"P-value: {p_seqkit:.4f}")

# Separate the meandepth data for each group
df_all = df[df['method'] == 'all']['meandepth']
df_bbnorm = df[df['method'] == 'bbnorm']['meandepth']
df_seqkit = df[df['method'] == 'seqkit']['meandepth']

print("="*40)
print("  Comparing meandepth for Each Method to 'all'")
print("="*40)

# --- Test 1: 'all' vs 'bbnorm' ---
stat_bbnorm, p_bbnorm = stats.ttest_ind(df_all, df_bbnorm)
print(f"Comparison: 'all' vs 'bbnorm'")
print(f"P-value: {p_bbnorm:.4f}")

# --- Test 2: 'all' vs 'seqkit' ---
stat_seqkit, p_seqkit = stats.ttest_ind(df_all, df_seqkit)
print(f"\nComparison: 'all' vs 'seqkit'")
print(f"P-value: {p_seqkit:.4f}")


# Separate the meanmapq data for each group
df_all = df[df['method'] == 'all']['meanbaseq']
df_bbnorm = df[df['method'] == 'bbnorm']['meanbaseq']
df_seqkit = df[df['method'] == 'seqkit']['meanbaseq']

print("="*40)
print("  Comparing meanbaseq for Each Method to 'all'")
print("="*40)

# --- Test 1: 'all' vs 'bbnorm' ---
stat_bbnorm, p_bbnorm = stats.ttest_ind(df_all, df_bbnorm)
print(f"Comparison: 'all' vs 'bbnorm'")
print(f"P-value: {p_bbnorm:.4f}")

# --- Test 2: 'all' vs 'seqkit' ---
stat_seqkit, p_seqkit = stats.ttest_ind(df_all, df_seqkit)
print(f"\nComparison: 'all' vs 'seqkit'")
print(f"P-value: {p_seqkit:.4f}")

# Separate the meanmapq data for each group
df_all = df[df['method'] == 'all']['meanmapq']
df_bbnorm = df[df['method'] == 'bbnorm']['meanmapq']
df_seqkit = df[df['method'] == 'seqkit']['meanmapq']

print("="*40)
print("  Comparing meanmapq for Each Method to 'all'")
print("="*40)

# --- Test 1: 'all' vs 'bbnorm' ---
stat_bbnorm, p_bbnorm = stats.ttest_ind(df_all, df_bbnorm)
print(f"Comparison: 'all' vs 'bbnorm'")
print(f"P-value: {p_bbnorm:.4f}")

# --- Test 2: 'all' vs 'seqkit' ---
stat_seqkit, p_seqkit = stats.ttest_ind(df_all, df_seqkit)
print(f"\nComparison: 'all' vs 'seqkit'")
print(f"P-value: {p_seqkit:.4f}")

# --- Bonferroni Correction ---
num_tests = 10
alpha = 0.05
bonferroni_alpha = alpha / num_tests

print("\n" + "="*40)
print("  Interpreting the Results")
print("="*40)
print(f"With a Bonferroni correction for {num_tests} tests, the adjusted significance level is {bonferroni_alpha}.")
print(f"Any p-value below {bonferroni_alpha} is considered statistically significant.")