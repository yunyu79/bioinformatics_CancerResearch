#%%
import pandas as pd
import numpy as np
# import os

#%% 
# Read all .tsv count files into Python as dataframe
dir= "/Users/yongxin/Documents/Research/HPC/deSeq2/Rworkstation/data/"

# ESCC, TCGA-IC-A6RF
tumor_00cf_1 = pd.read_csv(dir + "00cf_tumor_telescope-telescope_report.tsv", header = None, sep = "\t")
normal_b855_1 = pd.read_csv(dir + "b855_normal_telescope-telescope_report.tsv", header = None, sep = "\t")

# EAC, TCGA-L5-A43C, Mucinous adenocarcinoma counts for EAC
normal_8ec0_2 = pd.read_csv(dir + "8ec0_normal_telescope-telescope_report.tsv", header = None, sep = "\t")
tumor_b9d0_2 = pd.read_csv(dir + "b9d0_tumor_telescope-telescope_report.tsv", header = None, sep = "\t")

# EAC @SM:TCGA-L5-A4OG
tumor_0166_3 = pd.read_csv(dir + "0166_tumor_telescope-telescope_report.tsv", header = None, sep = "\t")
normal_ba94_3 = pd.read_csv(dir + "ba94_normal_telescope-telescope_report.tsv", header = None, sep = "\t")

# ESCC, SM:TCGA-L5-A4OM Basaloid squamous cell carcinoma
normal_9281_4 = pd.read_csv(dir + "9281_normal_telescope-telescope_report.tsv", header = None, sep = "\t")
tumor_cf39_4 = pd.read_csv(dir + "cf39_tumor_telescope-telescope_report.tsv", header = None, sep = "\t")

# EAC, TCGA-L5-A4OQ
normal_9cb3_5 = pd.read_csv(dir + "9cb3_normal_telescope-telescope_report.tsv", header = None, sep = "\t")
tumor_1421_5 = pd.read_csv(dir + "1421_tumor_telescope-telescope_report.tsv", header = None, sep = "\t")

# EAC, TCGA-V5-A7RE
normal_7a90_6 = pd.read_csv(dir + "7a90_normal_telescope-telescope_report.tsv", header = None, sep = "\t")
tumor_4225_6 = pd.read_csv(dir + "4225_tumor_telescope-telescope_report.tsv", header = None, sep = "\t")

# EAC, TCGA-IC-A6RE
tumor_6c97_7 = pd.read_csv(dir + "6c97_tumor_telescope-telescope_report.tsv", header = None, sep = "\t")
normal_e84e_7 = pd.read_csv(dir + "e84e_normal_telescope-telescope_report.tsv", header = None, sep = "\t")

# ESCC, TCGA-IG-A3I8
tumor_40ac_8 = pd.read_csv(dir + "40ac_tumor_telescope-telescope_report.tsv", header = None, sep = "\t")
normal_a3ba_8 = pd.read_csv(dir + "a3ba_normal_telescope-telescope_report.tsv", header = None, sep = "\t")

# EAC, TCGA-L5-A4OF
tumor_43e6_9 = pd.read_csv(dir + "43e6_tumor_telescope-telescope_report.tsv", header = None, sep = "\t")
normal_3562_9 = pd.read_csv(dir + "3562_normal_telescope-telescope_report.tsv", header = None, sep = "\t")

# EAC, TCGA-L5-A4OJ
tumor_5acb_10 = pd.read_csv(dir + "5acb_tumor_telescope-telescope_report.tsv", header = None, sep = "\t")
normal_caf4_10 = pd.read_csv(dir + "caf4_normal_telescope-telescope_report.tsv", header = None, sep = "\t")

# EAC, TCGA-L5-A4OO
normal_63d3_11 = pd.read_csv(dir + "63d3_normal_telescope-telescope_report.tsv", header = None, sep = "\t")
tumor_6259_11 = pd.read_csv(dir + "6259_tumor_telescope-telescope_report.tsv", header = None, sep = "\t")

# EAC, TCGA-L5-A4OR
normal_6d67_12 = pd.read_csv(dir + "6d67_normal_telescope-telescope_report.tsv", header = None, sep = "\t")
tumor_bdb3_12 = pd.read_csv(dir + "bdb3_tumor_telescope-telescope_report.tsv", header = None, sep = "\t")

# EAC, TCGA-V5-AASX
tumor_2f8e_13 = pd.read_csv(dir + "2f8e_tumor_telescope-telescope_report.tsv", header = None, sep = "\t")
normal_62b6_13 = pd.read_csv(dir + "62b6_normal_telescope-telescope_report.tsv", header = None, sep = "\t")

#%% 
# transform all tumor samples in a loop with list to fit in DeSeq2 countdata matrix standards
tumor_list = [tumor_00cf_1, tumor_b9d0_2, tumor_0166_3, tumor_cf39_4, tumor_1421_5, tumor_4225_6, tumor_6c97_7, 
              tumor_40ac_8, tumor_43e6_9, tumor_5acb_10, tumor_6259_11, tumor_bdb3_12, tumor_2f8e_13]

for sample_tumor in tumor_list:
    sample_tumor.drop([0, 2], inplace= True) # drop the irrelerant rows: 1st and 'no_feature'
    sample_tumor.columns = sample_tumor.iloc[0] # set the new first row as header
    sample_tumor.drop(index=sample_tumor.index[0], inplace=True)
    # tumor_sample = tumor_sample.iloc[1:] # keep all rows except the first row (new header)
    # tumor_trans_list.append(sample_tumor.iloc[:,[0,2]]) # keep the final_count column
    sample_tumor.drop(['transcript_length', 'final_conf', 'final_prop', 'init_aligned', 'unique_count', 'init_best', 'init_best_random', 'init_best_avg', 'init_prop', np.nan], axis=1, inplace = True)
    print(sample_tumor.head(3))

# transform all normal samples as dataframes to fit in DeSeq2 countdata matrix
normal_list = [normal_b855_1, normal_8ec0_2, normal_ba94_3, normal_9281_4, normal_9cb3_5, normal_7a90_6, normal_e84e_7,
               normal_a3ba_8, normal_3562_9, normal_caf4_10, normal_63d3_11, normal_6d67_12, normal_62b6_13]

for sample_normal in normal_list:
    sample_normal.drop([0, 2], inplace= True) # drop the first row
    sample_normal.columns = sample_normal.iloc[0] # set the new first row as header
    sample_normal.drop(index=sample_normal.index[0], inplace=True)
    # tumor_sample = tumor_sample.iloc[1:] # keep all rows except the first row (new header)
    # tumor_trans_list.append(sample_tumor.iloc[:,[0,2]]) # keep the final_count column
    sample_normal.drop(['transcript_length', 'final_conf', 'final_prop', 'init_aligned', 'unique_count', 'init_best', 'init_best_random', 'init_best_avg', 'init_prop', np.nan], axis=1, inplace = True)
    print(sample_normal.head(3))


# %%
# merge two dataframes into one with same row_name
samples_list = tumor_list + normal_list

# 
merged_df = samples_list[0]

for index, sample2merge in enumerate(samples_list[1:]):
    suffix = f"_df{index+1}"
    merged_df = merged_df.merge(sample2merge, how = "outer", on = "transcript", suffixes=('', suffix))

merged_df

merged_df.columns = ["transcript", "tumor_00cf_1", "tumor_b9d0_2", "tumor_0166_3", "tumor_cf39_4", "tumor_1421_5", "tumor_4225_6", 
                     "tumor_6c97_7", "tumor_40ac_8", "tumor_43e6_9", "tumor_5acb_10", "tumor_6259_11", "tumor_bdb3_12", "tumor_2f8e_13", 
                     "normal_b855_1", "normal_8ec0_2", "normal_ba94_3", "normal_9281_4", "normal_9cb3_5", "normal_7a90_6", 
                     "normal_e84e_7", "normal_a3ba_8", "normal_3562_9", "normal_caf4_10", "normal_63d3_11", "normal_6d67_12", "normal_62b6_13"]
merged_df.head(4)

# Fill NA with 0 to convert datatype from string to int
merged_df.fillna(0, inplace=True)
df_types = merged_df.applymap(type)

for column_name in merged_df.columns[1:]:
    merged_df[column_name] = merged_df[column_name].astype('int64')

df_types = merged_df.applymap(type)

# Replace 0 with NA to drop rows with all N/A on count-columns
merged_df.replace(0, np.nan, inplace=True)
merged_df.dropna(axis= 0, how = 'all', subset = merged_df.columns[1:], inplace= True)

# N/A --> 0 & output
merged_df.fillna(0, inplace = True)
merged_df.to_csv("finalCounts_pool.csv", index = False)

# filtered_df = merged_df[merged_df.iloc[:,1:].sum(axis=1)>10]
# filtered_df.to_csv("finalCounts_pool_filted.csv", index = False)

# %%
