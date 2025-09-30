import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# file header
# sample,bbnorm_bbnorm_10_cpu,bbnorm_bbnorm_10_time,bbnorm_bbnorm_10_memory,bbnorm_fastp_10_cpu,bbnorm_fastp_10_time,bbnorm_fastp_10_memory,bbnorm_fastp_30_cpu,bbnorm_fastp_30_time,bbnorm_fastp_30_memory,bbnorm_fastp_50_cpu,bbnorm_fastp_50_time,bbnorm_fastp_50_memory,bbnorm_fastp_100_cpu,bbnorm_fastp_100_time,bbnorm_fastp_100_memory,bbnorm_fastp_200_cpu,bbnorm_fastp_200_time,bbnorm_fastp_200_memory,bbnorm_fastp_300_cpu,bbnorm_fastp_300_time,bbnorm_fastp_300_memory,bbnorm_fastp_400_cpu,bbnorm_fastp_400_time,bbnorm_fastp_400_memory,bbnorm_fastp_500_cpu,bbnorm_fastp_500_time,bbnorm_fastp_500_memory,bbnorm_fastp_1000_cpu,bbnorm_fastp_1000_time,bbnorm_fastp_1000_memory,bbnorm_fastp_5000_cpu,bbnorm_fastp_5000_time,bbnorm_fastp_5000_memory,bbnorm_fastp_10000_cpu,bbnorm_fastp_10000_time,bbnorm_fastp_10000_memory,bbnorm_bbmap_10_cpu,bbnorm_bbmap_10_time,bbnorm_bbmap_10_memory,bbnorm_bbmap_30_cpu,bbnorm_bbmap_30_time,bbnorm_bbmap_30_memory,bbnorm_bbmap_50_cpu,bbnorm_bbmap_50_time,bbnorm_bbmap_50_memory,bbnorm_bbmap_100_cpu,bbnorm_bbmap_100_time,bbnorm_bbmap_100_memory,bbnorm_bbmap_200_cpu,bbnorm_bbmap_200_time,bbnorm_bbmap_200_memory,bbnorm_bbmap_300_cpu,bbnorm_bbmap_300_time,bbnorm_bbmap_300_memory,bbnorm_bbmap_400_cpu,bbnorm_bbmap_400_time,bbnorm_bbmap_400_memory,bbnorm_bbmap_500_cpu,bbnorm_bbmap_500_time,bbnorm_bbmap_500_memory,bbnorm_bbmap_1000_cpu,bbnorm_bbmap_1000_time,bbnorm_bbmap_1000_memory,bbnorm_bbmap_5000_cpu,bbnorm_bbmap_5000_time,bbnorm_bbmap_5000_memory,bbnorm_bbmap_10000_cpu,bbnorm_bbmap_10000_time,bbnorm_bbmap_10000_memory,bbnorm_ivar_consensus_consensus_10_cpu,bbnorm_ivar_consensus_consensus_10_time,bbnorm_ivar_consensus_consensus_10_memory,bbnorm_ivar_consensus_consensus_30_cpu,bbnorm_ivar_consensus_consensus_30_time,bbnorm_ivar_consensus_consensus_30_memory,bbnorm_ivar_consensus_consensus_50_cpu,bbnorm_ivar_consensus_consensus_50_time,bbnorm_ivar_consensus_consensus_50_memory,bbnorm_ivar_consensus_consensus_100_cpu,bbnorm_ivar_consensus_consensus_100_time,bbnorm_ivar_consensus_consensus_100_memory,bbnorm_ivar_consensus_consensus_200_cpu,bbnorm_ivar_consensus_consensus_200_time,bbnorm_ivar_consensus_consensus_200_memory,bbnorm_ivar_consensus_consensus_300_cpu,bbnorm_ivar_consensus_consensus_300_time,bbnorm_ivar_consensus_consensus_300_memory,bbnorm_ivar_consensus_consensus_400_cpu,bbnorm_ivar_consensus_consensus_400_time,bbnorm_ivar_consensus_consensus_400_memory,bbnorm_ivar_consensus_consensus_500_cpu,bbnorm_ivar_consensus_consensus_500_time,bbnorm_ivar_consensus_consensus_500_memory,bbnorm_ivar_consensus_consensus_1000_cpu,bbnorm_ivar_consensus_consensus_1000_time,bbnorm_ivar_consensus_consensus_1000_memory,bbnorm_ivar_consensus_consensus_5000_cpu,bbnorm_ivar_consensus_consensus_5000_time,bbnorm_ivar_consensus_consensus_5000_memory,bbnorm_ivar_consensus_consensus_10000_cpu,bbnorm_ivar_consensus_consensus_10000_time,bbnorm_ivar_consensus_consensus_10000_memory,bbnorm_ivar_consensus_mpileup_10_cpu,bbnorm_ivar_consensus_mpileup_10_time,bbnorm_ivar_consensus_mpileup_10_memory,bbnorm_ivar_consensus_mpileup_30_cpu,bbnorm_ivar_consensus_mpileup_30_time,bbnorm_ivar_consensus_mpileup_30_memory,bbnorm_ivar_consensus_mpileup_50_cpu,bbnorm_ivar_consensus_mpileup_50_time,bbnorm_ivar_consensus_mpileup_50_memory,bbnorm_ivar_consensus_mpileup_100_cpu,bbnorm_ivar_consensus_mpileup_100_time,bbnorm_ivar_consensus_mpileup_100_memory,bbnorm_ivar_consensus_mpileup_200_cpu,bbnorm_ivar_consensus_mpileup_200_time,bbnorm_ivar_consensus_mpileup_200_memory,bbnorm_ivar_consensus_mpileup_300_cpu,bbnorm_ivar_consensus_mpileup_300_time,bbnorm_ivar_consensus_mpileup_300_memory,bbnorm_ivar_consensus_mpileup_400_cpu,bbnorm_ivar_consensus_mpileup_400_time,bbnorm_ivar_consensus_mpileup_400_memory,bbnorm_ivar_consensus_mpileup_500_cpu,bbnorm_ivar_consensus_mpileup_500_time,bbnorm_ivar_consensus_mpileup_500_memory,bbnorm_ivar_consensus_mpileup_1000_cpu,bbnorm_ivar_consensus_mpileup_1000_time,bbnorm_ivar_consensus_mpileup_1000_memory,bbnorm_ivar_consensus_mpileup_5000_cpu,bbnorm_ivar_consensus_mpileup_5000_time,bbnorm_ivar_consensus_mpileup_5000_memory,bbnorm_ivar_consensus_mpileup_10000_cpu,bbnorm_ivar_consensus_mpileup_10000_time,bbnorm_ivar_consensus_mpileup_10000_memory,bbnorm_ivar_trim_10_cpu,bbnorm_ivar_trim_10_time,bbnorm_ivar_trim_10_memory,bbnorm_ivar_trim_30_cpu,bbnorm_ivar_trim_30_time,bbnorm_ivar_trim_30_memory,bbnorm_ivar_trim_50_cpu,bbnorm_ivar_trim_50_time,bbnorm_ivar_trim_50_memory,bbnorm_ivar_trim_100_cpu,bbnorm_ivar_trim_100_time,bbnorm_ivar_trim_100_memory,bbnorm_ivar_trim_200_cpu,bbnorm_ivar_trim_200_time,bbnorm_ivar_trim_200_memory,bbnorm_ivar_trim_300_cpu,bbnorm_ivar_trim_300_time,bbnorm_ivar_trim_300_memory,bbnorm_ivar_trim_400_cpu,bbnorm_ivar_trim_400_time,bbnorm_ivar_trim_400_memory,bbnorm_ivar_trim_500_cpu,bbnorm_ivar_trim_500_time,bbnorm_ivar_trim_500_memory,bbnorm_ivar_trim_1000_cpu,bbnorm_ivar_trim_1000_time,bbnorm_ivar_trim_1000_memory,bbnorm_ivar_trim_5000_cpu,bbnorm_ivar_trim_5000_time,bbnorm_ivar_trim_5000_memory,bbnorm_ivar_trim_10000_cpu,bbnorm_ivar_trim_10000_time,bbnorm_ivar_trim_10000_memory,bbnorm_sort_10_cpu,bbnorm_sort_10_time,bbnorm_sort_10_memory,bbnorm_sort_30_cpu,bbnorm_sort_30_time,bbnorm_sort_30_memory,bbnorm_sort_50_cpu,bbnorm_sort_50_time,bbnorm_sort_50_memory,bbnorm_sort_100_cpu,bbnorm_sort_100_time,bbnorm_sort_100_memory,bbnorm_sort_200_cpu,bbnorm_sort_200_time,bbnorm_sort_200_memory,bbnorm_sort_300_cpu,bbnorm_sort_300_time,bbnorm_sort_300_memory,bbnorm_sort_400_cpu,bbnorm_sort_400_time,bbnorm_sort_400_memory,bbnorm_sort_500_cpu,bbnorm_sort_500_time,bbnorm_sort_500_memory,bbnorm_sort_1000_cpu,bbnorm_sort_1000_time,bbnorm_sort_1000_memory,bbnorm_sort_5000_cpu,bbnorm_sort_5000_time,bbnorm_sort_5000_memory,bbnorm_sort_10000_cpu,bbnorm_sort_10000_time,bbnorm_sort_10000_memory,bbnorm_sort2_10_cpu,bbnorm_sort2_10_time,bbnorm_sort2_10_memory,bbnorm_sort2_30_cpu,bbnorm_sort2_30_time,bbnorm_sort2_30_memory,bbnorm_sort2_50_cpu,bbnorm_sort2_50_time,bbnorm_sort2_50_memory,bbnorm_sort2_100_cpu,bbnorm_sort2_100_time,bbnorm_sort2_100_memory,bbnorm_sort2_200_cpu,bbnorm_sort2_200_time,bbnorm_sort2_200_memory,bbnorm_sort2_300_cpu,bbnorm_sort2_300_time,bbnorm_sort2_300_memory,bbnorm_sort2_400_cpu,bbnorm_sort2_400_time,bbnorm_sort2_400_memory,bbnorm_sort2_500_cpu,bbnorm_sort2_500_time,bbnorm_sort2_500_memory,bbnorm_sort2_1000_cpu,bbnorm_sort2_1000_time,bbnorm_sort2_1000_memory,bbnorm_sort2_5000_cpu,bbnorm_sort2_5000_time,bbnorm_sort2_5000_memory,bbnorm_sort2_10000_cpu,bbnorm_sort2_10000_time,bbnorm_sort2_10000_memory,seqkit_seqkit_990_R1_cpu,seqkit_seqkit_990_R1_time,seqkit_seqkit_990_R1_memory,seqkit_seqkit_990_R2_cpu,seqkit_seqkit_990_R2_time,seqkit_seqkit_990_R2_memory,seqkit_seqkit_2970_R1_cpu,seqkit_seqkit_2970_R1_time,seqkit_seqkit_2970_R1_memory,seqkit_seqkit_2970_R2_cpu,seqkit_seqkit_2970_R2_time,seqkit_seqkit_2970_R2_memory,seqkit_seqkit_4951_R1_cpu,seqkit_seqkit_4951_R1_time,seqkit_seqkit_4951_R1_memory,seqkit_seqkit_4951_R2_cpu,seqkit_seqkit_4951_R2_time,seqkit_seqkit_4951_R2_memory,seqkit_seqkit_9902_R1_cpu,seqkit_seqkit_9902_R1_time,seqkit_seqkit_9902_R1_memory,seqkit_seqkit_9902_R2_cpu,seqkit_seqkit_9902_R2_time,seqkit_seqkit_9902_R2_memory,seqkit_seqkit_14852_R1_cpu,seqkit_seqkit_14852_R1_time,seqkit_seqkit_14852_R1_memory,seqkit_seqkit_14852_R2_cpu,seqkit_seqkit_14852_R2_time,seqkit_seqkit_14852_R2_memory,seqkit_seqkit_19803_R1_cpu,seqkit_seqkit_19803_R1_time,seqkit_seqkit_19803_R1_memory,seqkit_seqkit_19803_R2_cpu,seqkit_seqkit_19803_R2_time,seqkit_seqkit_19803_R2_memory,seqkit_seqkit_29705_R1_cpu,seqkit_seqkit_29705_R1_time,seqkit_seqkit_29705_R1_memory,seqkit_seqkit_29705_R2_cpu,seqkit_seqkit_29705_R2_time,seqkit_seqkit_29705_R2_memory,seqkit_seqkit_49508_R1_cpu,seqkit_seqkit_49508_R1_time,seqkit_seqkit_49508_R1_memory,seqkit_seqkit_49508_R2_cpu,seqkit_seqkit_49508_R2_time,seqkit_seqkit_49508_R2_memory,seqkit_seqkit_99017_R1_cpu,seqkit_seqkit_99017_R1_time,seqkit_seqkit_99017_R1_memory,seqkit_seqkit_99017_R2_cpu,seqkit_seqkit_99017_R2_time,seqkit_seqkit_99017_R2_memory,seqkit_seqkit_495083_R1_cpu,seqkit_seqkit_495083_R1_time,seqkit_seqkit_495083_R1_memory,seqkit_seqkit_495083_R2_cpu,seqkit_seqkit_495083_R2_time,seqkit_seqkit_495083_R2_memory,seqkit_seqkit_990166_R1_cpu,seqkit_seqkit_990166_R1_time,seqkit_seqkit_990166_R1_memory,seqkit_seqkit_990166_R2_cpu,seqkit_seqkit_990166_R2_time,seqkit_seqkit_990166_R2_memory,seqkit_fastp_10_cpu,seqkit_fastp_10_time,seqkit_fastp_10_memory,seqkit_fastp_30_cpu,seqkit_fastp_30_time,seqkit_fastp_30_memory,seqkit_fastp_50_cpu,seqkit_fastp_50_time,seqkit_fastp_50_memory,seqkit_fastp_100_cpu,seqkit_fastp_100_time,seqkit_fastp_100_memory,seqkit_fastp_200_cpu,seqkit_fastp_200_time,seqkit_fastp_200_memory,seqkit_fastp_300_cpu,seqkit_fastp_300_time,seqkit_fastp_300_memory,seqkit_fastp_400_cpu,seqkit_fastp_400_time,seqkit_fastp_400_memory,seqkit_fastp_500_cpu,seqkit_fastp_500_time,seqkit_fastp_500_memory,seqkit_fastp_1000_cpu,seqkit_fastp_1000_time,seqkit_fastp_1000_memory,seqkit_fastp_5000_cpu,seqkit_fastp_5000_time,seqkit_fastp_5000_memory,seqkit_fastp_10000_cpu,seqkit_fastp_10000_time,seqkit_fastp_10000_memory,seqkit_bbmap_10_cpu,seqkit_bbmap_10_time,seqkit_bbmap_10_memory,seqkit_bbmap_30_cpu,seqkit_bbmap_30_time,seqkit_bbmap_30_memory,seqkit_bbmap_50_cpu,seqkit_bbmap_50_time,seqkit_bbmap_50_memory,seqkit_bbmap_100_cpu,seqkit_bbmap_100_time,seqkit_bbmap_100_memory,seqkit_bbmap_200_cpu,seqkit_bbmap_200_time,seqkit_bbmap_200_memory,seqkit_bbmap_300_cpu,seqkit_bbmap_300_time,seqkit_bbmap_300_memory,seqkit_bbmap_400_cpu,seqkit_bbmap_400_time,seqkit_bbmap_400_memory,seqkit_bbmap_500_cpu,seqkit_bbmap_500_time,seqkit_bbmap_500_memory,seqkit_bbmap_1000_cpu,seqkit_bbmap_1000_time,seqkit_bbmap_1000_memory,seqkit_bbmap_5000_cpu,seqkit_bbmap_5000_time,seqkit_bbmap_5000_memory,seqkit_bbmap_10000_cpu,seqkit_bbmap_10000_time,seqkit_bbmap_10000_memory,seqkit_ivar_consensus_consensus_10_cpu,seqkit_ivar_consensus_consensus_10_time,seqkit_ivar_consensus_consensus_10_memory,seqkit_ivar_consensus_consensus_30_cpu,seqkit_ivar_consensus_consensus_30_time,seqkit_ivar_consensus_consensus_30_memory,seqkit_ivar_consensus_consensus_50_cpu,seqkit_ivar_consensus_consensus_50_time,seqkit_ivar_consensus_consensus_50_memory,seqkit_ivar_consensus_consensus_100_cpu,seqkit_ivar_consensus_consensus_100_time,seqkit_ivar_consensus_consensus_100_memory,seqkit_ivar_consensus_consensus_200_cpu,seqkit_ivar_consensus_consensus_200_time,seqkit_ivar_consensus_consensus_200_memory,seqkit_ivar_consensus_consensus_300_cpu,seqkit_ivar_consensus_consensus_300_time,seqkit_ivar_consensus_consensus_300_memory,seqkit_ivar_consensus_consensus_400_cpu,seqkit_ivar_consensus_consensus_400_time,seqkit_ivar_consensus_consensus_400_memory,seqkit_ivar_consensus_consensus_500_cpu,seqkit_ivar_consensus_consensus_500_time,seqkit_ivar_consensus_consensus_500_memory,seqkit_ivar_consensus_consensus_1000_cpu,seqkit_ivar_consensus_consensus_1000_time,seqkit_ivar_consensus_consensus_1000_memory,seqkit_ivar_consensus_consensus_5000_cpu,seqkit_ivar_consensus_consensus_5000_time,seqkit_ivar_consensus_consensus_5000_memory,seqkit_ivar_consensus_consensus_10000_cpu,seqkit_ivar_consensus_consensus_10000_time,seqkit_ivar_consensus_consensus_10000_memory,seqkit_ivar_consensus_mpileup_10_cpu,seqkit_ivar_consensus_mpileup_10_time,seqkit_ivar_consensus_mpileup_10_memory,seqkit_ivar_consensus_mpileup_30_cpu,seqkit_ivar_consensus_mpileup_30_time,seqkit_ivar_consensus_mpileup_30_memory,seqkit_ivar_consensus_mpileup_50_cpu,seqkit_ivar_consensus_mpileup_50_time,seqkit_ivar_consensus_mpileup_50_memory,seqkit_ivar_consensus_mpileup_100_cpu,seqkit_ivar_consensus_mpileup_100_time,seqkit_ivar_consensus_mpileup_100_memory,seqkit_ivar_consensus_mpileup_200_cpu,seqkit_ivar_consensus_mpileup_200_time,seqkit_ivar_consensus_mpileup_200_memory,seqkit_ivar_consensus_mpileup_300_cpu,seqkit_ivar_consensus_mpileup_300_time,seqkit_ivar_consensus_mpileup_300_memory,seqkit_ivar_consensus_mpileup_400_cpu,seqkit_ivar_consensus_mpileup_400_time,seqkit_ivar_consensus_mpileup_400_memory,seqkit_ivar_consensus_mpileup_500_cpu,seqkit_ivar_consensus_mpileup_500_time,seqkit_ivar_consensus_mpileup_500_memory,seqkit_ivar_consensus_mpileup_1000_cpu,seqkit_ivar_consensus_mpileup_1000_time,seqkit_ivar_consensus_mpileup_1000_memory,seqkit_ivar_consensus_mpileup_5000_cpu,seqkit_ivar_consensus_mpileup_5000_time,seqkit_ivar_consensus_mpileup_5000_memory,seqkit_ivar_consensus_mpileup_10000_cpu,seqkit_ivar_consensus_mpileup_10000_time,seqkit_ivar_consensus_mpileup_10000_memory,seqkit_ivar_trim_10_cpu,seqkit_ivar_trim_10_time,seqkit_ivar_trim_10_memory,seqkit_ivar_trim_30_cpu,seqkit_ivar_trim_30_time,seqkit_ivar_trim_30_memory,seqkit_ivar_trim_50_cpu,seqkit_ivar_trim_50_time,seqkit_ivar_trim_50_memory,seqkit_ivar_trim_100_cpu,seqkit_ivar_trim_100_time,seqkit_ivar_trim_100_memory,seqkit_ivar_trim_200_cpu,seqkit_ivar_trim_200_time,seqkit_ivar_trim_200_memory,seqkit_ivar_trim_300_cpu,seqkit_ivar_trim_300_time,seqkit_ivar_trim_300_memory,seqkit_ivar_trim_400_cpu,seqkit_ivar_trim_400_time,seqkit_ivar_trim_400_memory,seqkit_ivar_trim_500_cpu,seqkit_ivar_trim_500_time,seqkit_ivar_trim_500_memory,seqkit_ivar_trim_1000_cpu,seqkit_ivar_trim_1000_time,seqkit_ivar_trim_1000_memory,seqkit_ivar_trim_5000_cpu,seqkit_ivar_trim_5000_time,seqkit_ivar_trim_5000_memory,seqkit_ivar_trim_10000_cpu,seqkit_ivar_trim_10000_time,seqkit_ivar_trim_10000_memory,seqkit_sort_10_cpu,seqkit_sort_10_time,seqkit_sort_10_memory,seqkit_sort_30_cpu,seqkit_sort_30_time,seqkit_sort_30_memory,seqkit_sort_50_cpu,seqkit_sort_50_time,seqkit_sort_50_memory,seqkit_sort_100_cpu,seqkit_sort_100_time,seqkit_sort_100_memory,seqkit_sort_200_cpu,seqkit_sort_200_time,seqkit_sort_200_memory,seqkit_sort_300_cpu,seqkit_sort_300_time,seqkit_sort_300_memory,seqkit_sort_400_cpu,seqkit_sort_400_time,seqkit_sort_400_memory,seqkit_sort_500_cpu,seqkit_sort_500_time,seqkit_sort_500_memory,seqkit_sort_1000_cpu,seqkit_sort_1000_time,seqkit_sort_1000_memory,seqkit_sort_5000_cpu,seqkit_sort_5000_time,seqkit_sort_5000_memory,seqkit_sort_10000_cpu,seqkit_sort_10000_time,seqkit_sort_10000_memory,seqkit_sort2_10_cpu,seqkit_sort2_10_time,seqkit_sort2_10_memory,seqkit_sort2_30_cpu,seqkit_sort2_30_time,seqkit_sort2_30_memory,seqkit_sort2_50_cpu,seqkit_sort2_50_time,seqkit_sort2_50_memory,seqkit_sort2_100_cpu,seqkit_sort2_100_time,seqkit_sort2_100_memory,seqkit_sort2_200_cpu,seqkit_sort2_200_time,seqkit_sort2_200_memory,seqkit_sort2_300_cpu,seqkit_sort2_300_time,seqkit_sort2_300_memory,seqkit_sort2_400_cpu,seqkit_sort2_400_time,seqkit_sort2_400_memory,seqkit_sort2_500_cpu,seqkit_sort2_500_time,seqkit_sort2_500_memory,seqkit_sort2_1000_cpu,seqkit_sort2_1000_time,seqkit_sort2_1000_memory,seqkit_sort2_5000_cpu,seqkit_sort2_5000_time,seqkit_sort2_5000_memory,seqkit_sort2_10000_cpu,seqkit_sort2_10000_time,seqkit_sort2_10000_memory,raw_fastp_all_cpu,raw_fastp_all_time,raw_fastp_all_memory,raw_bbmap_all_cpu,raw_bbmap_all_time,raw_bbmap_all_memory,raw_ivar_consensus_consensus_all_cpu,raw_ivar_consensus_consensus_all_time,raw_ivar_consensus_consensus_all_memory,raw_ivar_consensus_mpileup_all_cpu,raw_ivar_consensus_mpileup_all_time,raw_ivar_consensus_mpileup_all_memory,raw_ivar_trim_all_cpu,raw_ivar_trim_all_time,raw_ivar_trim_all_memory,raw_sort_all_cpu,raw_sort_all_time,raw_sort_all_memory,raw_sort2_all_cpu,raw_sort2_all_time,raw_sort2_all_memory

def parse_col_name(col_name):
    """Parses a column name into its constituent parts."""
    parts = col_name.split('_')

    method = ""
    if "bbnorm" in col_name:
        method = "bbnorm"
    elif "seqkit" in col_name:
        method = "seqkit"
    elif "raw" in col_name:
        method = "raw"
    else:
        print(parts)
        print(col_name)
        exit()

    step = ""
    for stp in ["fastp", "bbmap", "ivar_trim", "sort2", "ivar_consensus_mpileup", "ivar_consensus_consensus"]:
        if stp in col_name:
            step = stp
            break

    if parts[1] == "sort":
        step = "sort1"

    depth = parts[-2]
    resource = parts[-1]

    if "bbnorm_bbnorm" in col_name:
        step = "bbnorm"
    elif "seqkit_seqkit" in col_name:
        step = f"seqkit_{parts[3]}"
        depth = parts[2]

    if not method:
        print(col_name)
        print(parts)
        exit()

    if not step:
        print(col_name)
        print(parts)
        exit()

    num_depth_mapping = {
        990166 : 10000,
        495083 : 5000,
        99017 : 1000,
        49508 : 500,
        29705 : 400,
        19803 : 300,
        14852 : 200,
        9902 : 100,
        4951 : 50,
        2970 : 30,
        990 : 10
        }
    
    if method == "seqkit" and "seqkit" in step :
        depth = num_depth_mapping.get(int(depth))

    return method, step, depth, resource

def time_to_seconds(time_str):
    """Converts a time string (H:M:S or M:S) to total seconds."""
    parts = str(time_str).split(':')
    if len(parts) == 3:
        return int(parts[0]) * 3600 + int(parts[1]) * 60 + float(parts[2])
    elif len(parts) == 2:
        return int(parts[0]) * 60 + float(parts[1])
    else:
        return float(parts[0])


df = pd.read_csv("data/total_resources.csv")

# 2. Reshape the DataFrame from wide to long format
df_long = df.melt(id_vars=['sample'], var_name='metric', value_name='value')

# 3. Parse the metric column into separate columns for easier filtering
parsed_metrics = df_long['metric'].apply(parse_col_name).apply(pd.Series)
parsed_metrics.columns = ['method', 'step', 'depth', 'resource']
df_long = pd.concat([df_long, parsed_metrics], axis=1)

steps = [
    "bbnorm",
    "seqkit_R1", 
    "seqkit_R2", 
    "fastp", 
    "bbmap", 
    "sort1", 
    "ivar_trim", 
    "sort2", 
    "ivar_consensus_mpileup", 
    "ivar_consensus_consensus"
    ]

methods = ['raw', 'bbnorm', 'seqkit']

boxplot_color_map = {
    'raw': '#0173B2',    # Blue
    'bbnorm': '#DE8F05', # Orange
    'seqkit': '#029E73'  # Green
}

desired_depths = ['10', '30', '50','100', '200', '300', '400', '500', '1000', '5000', '10000', 'all']
depth_order    = ['10', '30', '50','100', '200', '300', '400', '500', '1000', '5000', '10000', 'all']

# for stp in steps:
#     filtered_df = df_long[df_long['step'] == stp]
#     print(f"Step: {stp}")
#     print(filtered_df.head(), "\n")

# for method in methods:
#     filtered_df = df_long[df_long['method'] == method]
#     print(f"Step: {stp}")
#     print(filtered_df.head(), "\n")

df_long = df_long.drop(columns=['metric'])
df_long['method_depth'] = df_long['method'] + '_' + df_long['depth'].astype(str)

print("getting together cpu plot for all steps")
df_cpu = df_long[df_long['resource'] == 'cpu'].copy()

# The 'depth' column has mixed types, convert to string for filtering
df_cpu['depth_str'] = df_cpu['depth'].astype(str)
df_cpu['step'] = df_cpu['step'].str.replace(r'^seqkit_.*$', 'seqkit', regex=True)

# Filter for the specific depth levels of interest
plot_df = df_cpu[df_cpu['depth_str'].isin(desired_depths)].copy()

plt.figure(figsize=(12, 7))
sns.boxplot(
    data=plot_df, 
    x='depth_str', 
    y='value', 
    hue='method', 
    order=depth_order,
    palette=boxplot_color_map
    )
plt.title('Average CPU Usage by Method and Depth', fontsize=16)
plt.xlabel('Method', fontsize=12)
plt.ylabel('Average CPU Usage (percent)', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.6)
plt.savefig('cpu_usage_by_method_and_depth.png', dpi=300)
plt.close()


for stp in steps + ['seqkit']:
    print(f"getting together cpu plot for {stp}")
    stp_plot_df = plot_df[plot_df['step'] == stp].copy()

    plt.figure(figsize=(12, 7))
    sns.boxplot(
        data=stp_plot_df, 
        x='depth_str', 
        y='value', 
        hue='method', 
        order=depth_order,
        palette=boxplot_color_map
        )
    plt.title(f'Average CPU Usage by Method and Depth for {stp}', fontsize=16)
    plt.xlabel('Method', fontsize=12)
    plt.ylabel('Average CPU Usage (percent)', fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.savefig(f"cpu_usage_by_method_and_depth_for_{stp}.png", dpi=300)
    plt.close()


print("getting together memory plot for all steps")
df_memory = df_long[df_long['resource'] == 'memory'].copy()

# The 'depth' column has mixed types, convert to string for filtering
df_memory['depth_str'] = df_memory['depth'].astype(str)
df_memory['step'] = df_memory['step'].str.replace(r'^seqkit_.*$', 'seqkit', regex=True)


# Filter for the specific depth levels of interest
plot_df = df_memory[df_memory['depth_str'].isin(desired_depths)].copy()

plt.figure(figsize=(12, 7))
sns.boxplot(
    data=plot_df, 
    x='depth_str', 
    y='value', 
    hue='method', 
    order=depth_order,
    palette=boxplot_color_map
    )
plt.yscale('log')
plt.title('Average Memory Usage by Method and Depth', fontsize=16)
plt.xlabel('Method', fontsize=12)
plt.ylabel('Average Memory Usage (kbytes)', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.6)
plt.savefig('memory_usage_by_method_and_depth.png', dpi=300)
plt.close()


for stp in steps + ['seqkit']:
    print(f"getting together memory plot for {stp}")
    stp_plot_df = plot_df[plot_df['step'] == stp].copy()

    plt.figure(figsize=(12, 7))
    sns.boxplot(
        data=stp_plot_df, 
        x='depth_str', 
        y='value', 
        hue='method', 
        order=depth_order,
        palette=boxplot_color_map
        )
    plt.yscale('log')
    plt.title(f'Average Memory Usage by Method and Depth for {stp}', fontsize=16)
    plt.xlabel('Method', fontsize=12)
    plt.ylabel('Average Memory Usage (kbytes)', fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.savefig(f"memory_usage_by_method_and_depth_for_{stp}.png", dpi=300)
    plt.close()

print("getting together time plot for all steps")
df_time = df_long[df_long['resource'] == 'time'].copy()

# The 'depth' column has mixed types, convert to string for filtering
df_time['depth_str'] = df_time['depth'].astype(str)
df_time['value_seconds'] = df_time['value'].apply(time_to_seconds)

# Filter for the specific depth levels of interest
plot_df = df_time[df_time['depth_str'].isin(desired_depths)].copy()

plt.figure(figsize=(12, 7))
sns.boxplot(
    data=plot_df, 
    x='depth_str', 
    y='value_seconds', 
    hue='method', 
    order=depth_order,
    palette=boxplot_color_map
    )
plt.yscale('log')
plt.title('Average Time Usage by Method and Depth', fontsize=16)
plt.xlabel('Method', fontsize=12)
plt.ylabel('Average Time Usage (seconds)', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.6)
plt.savefig('time_usage_by_method_and_depth.png', dpi=300)
plt.close()


for stp in steps + ['seqkit']:
    print(f"getting together time plot for {stp}")
    stp_plot_df = plot_df[plot_df['step'] == stp].copy()

    plt.figure(figsize=(12, 7))
    sns.boxplot(
        data=stp_plot_df, 
        x='depth_str', 
        y='value_seconds', 
        hue='method', 
        order=depth_order,
        palette=boxplot_color_map
        )
    plt.yscale('log')
    plt.title(f'Average Time Length by Method and Depth for {stp}', fontsize=16)
    plt.xlabel('Method', fontsize=12)
    plt.ylabel('Average Time Usage (seconds)', fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.savefig(f"time_usage_by_method_and_depth_for_{stp}.png", dpi=300)
    plt.close()

step_order = [
    "bbnorm", 
    "seqkit_R1",
    "seqkit_R2",
    "fastp", 
    "bbmap", 
    "sort1", 
    "ivar_trim", 
    "sort2",
    "ivar_consensus_mpileup", 
    "ivar_consensus_consensus"
]

method_order = ['raw', 'bbnorm', 'seqkit']

method_depth_order = [f"{m}_{d}" for m in method_order for d in depth_order if f"{m}_{d}" in df_time['method_depth'].unique()]

df_time['method_depth'] = pd.Categorical(df_time['method_depth'], categories=method_depth_order, ordered=True)

avg_times = df_time.groupby(['method_depth', 'step'])['value_seconds'].mean().reset_index()
pivot_df = avg_times.pivot(index='method_depth', columns='step', values='value_seconds').fillna(0)
pivot_df = pivot_df[step_order]  # ensure correct order of columns

# Plot stacked bar chart
fig, ax = plt.subplots(figsize=(14, 7))
bottom = pd.Series([0]*len(pivot_df), index=pivot_df.index)

for step in step_order:
    ax.bar(pivot_df.index, pivot_df[step], bottom=bottom, label=step)
    bottom += pivot_df[step]

ax.set_ylabel("Average Time (seconds)")
ax.set_xlabel("Method_Depth")
ax.set_title("Average Time per Step by Method/Depth (seconds)")
ax.legend(title="Step", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.xticks(rotation=90, ha='right')
plt.tight_layout()

plt.savefig("time_usage_by_method_and_depth_bar_chart.png", dpi=300)
plt.close()

# Plot stacked bar chart
fig, ax = plt.subplots(figsize=(14, 7))
bottom = pd.Series([0]*len(pivot_df), index=pivot_df.index)

for step in step_order:
    ax.bar(pivot_df.index, pivot_df[step], bottom=bottom, label=step)
    bottom += pivot_df[step]

ax.set_ylabel("Average Time (seconds)")
ax.set_xlabel("Method_Depth")
ax.set_title("Average Time per Step by Method/Depth (seconds)")
ax.legend(title="Step", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.yscale('log')
plt.xticks(rotation=90, ha='right')
plt.tight_layout()

plt.savefig("time_usage_by_method_and_depth_bar_chart_logscale.png", dpi=300)
plt.close()