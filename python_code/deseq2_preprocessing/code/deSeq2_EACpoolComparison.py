#%%
import pandas as pd
import numpy as np
import os

#%%

dir_EAC_allTumor = "/Users/yongxin/Documents/Research/HPC/deSeq2/Rworkstation/6_telescope_EAC_tumor"
dir_EAC_allNormal = "/Users/yongxin/Documents/Research/HPC/deSeq2/Rworkstation/6_telescope_EAC_normal"
dataframes = {} # create a dictionary to store each sample's result as a dataframe (an element)

for filename in os.listdir(dir_EAC_allTumor):
    if filename.endswith(".tsv"):
        prefix = os.path.splitext(filename)[0]
        filepath = os.path.join(dir_EAC_allTumor,filename)
        
        df = pd.read_csv(filepath, header = None, sep = "\t",  encoding="utf-8")

        df.drop([0, 2], inplace= True) # drop the irrelerant rows: 1st and 'no_feature'
        df.columns = df.iloc[0] # set the new first row as header
        df.drop(index=df.index[0], inplace=True)
        # tumor_sample = tumor_sample.iloc[1:] # keep all rows except the first row (new header)
        # tumor_trans_list.append(sample_tumor.iloc[:,[0,2]]) # keep the final_count column
        df.drop(['transcript_length', 'final_conf', 'final_prop', 'init_aligned', 'unique_count', 'init_best', 'init_best_random', 'init_best_avg', 'init_prop', np.nan], axis=1, inplace = True)
        df.columns = ['transcript', "tumor_"+prefix]
        print(df.head(3))

        dataframes["tumor_"+prefix] = df  # add a new element to dictionary

print(dataframes.keys())
print(dataframes['tumor_1421_telescope'])

for filename in os.listdir(dir_EAC_allNormal):
    if filename.endswith(".tsv"):
        prefix = os.path.splitext(filename)[0]
        filepath = os.path.join(dir_EAC_allNormal,filename)
        
        df = pd.read_csv(filepath, header = None, sep = "\t",  encoding="utf-8")

        df.drop([0, 2], inplace= True) # drop the irrelerant rows: 1st and 'no_feature'
        df.columns = df.iloc[0] # set the new first row as header
        df.drop(index=df.index[0], inplace=True)
        # tumor_sample = tumor_sample.iloc[1:] # keep all rows except the first row (new header)
        # tumor_trans_list.append(sample_tumor.iloc[:,[0,2]]) # keep the final_count column
        df.drop(['transcript_length', 'final_conf', 'final_prop', 'init_aligned', 'unique_count', 'init_best', 'init_best_random', 'init_best_avg', 'init_prop', np.nan], axis=1, inplace = True)
        df.columns = ['transcript', "normal_"+prefix]
        print(df.head(3))

        dataframes["normal_"+prefix] = df 

print(dataframes.keys())

#%%
if dataframes:
    merged_df = next(iter(dataframes.values()))
    firstKey, firstValue = next(iter(dataframes.items()))
    dataframes.pop(firstKey)
    for name, df in dataframes.items():
        merged_df = merged_df.merge(df, how = "outer", on = "transcript")


# %%
# merge two dataframes into one with same row_name

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
merged_df.to_csv("countData_EAC_pool.csv", index = False)

# filtered_df = merged_df[merged_df.iloc[:,1:].sum(axis=1)>10]
# filtered_df.to_csv("finalCounts_pool_filted.csv", index = False)

# %%
