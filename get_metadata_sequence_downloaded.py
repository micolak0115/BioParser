# %%
import os
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
# %%
path_sra_metadata = "/BiO/Research/GeroPathway/Resources/RawData/Transcriptome/Killifish/GRZ/sra_metadata_run_info_GRZ.txt"
path_dir_sra_seq = "/BiO/Research/GeroPathway/Resources/RawData/Transcriptome/Killifish/GRZ"

# %%
df_sra_metadata = pd.read_csv(path_sra_metadata, sep="\t")
list_dir_sra_seq = os.listdir(path_dir_sra_seq)
list_dir_sra_seq = list(filter(lambda x: str(x).split(".")[-1] not in ["log", "txt"], list_dir_sra_seq))
df_sra_metdata_set_idx = df_sra_metadata.set_index("Run")
df_sra_metadata_filt = df_sra_metdata_set_idx.loc[list_dir_sra_seq]
list_sex = df_sra_metadata_filt["LibraryName"].to_list()
dict_sex_count = dict(Counter(list_sex))
for col in df_sra_metadata_filt.columns:
    list_sex = df_sra_metadata_filt[col].to_list()
    if len(set(list_sex)) < 20:
        dict_sex_count = dict(Counter(list_sex))
        plt.figure(figsize=(20,5))
        plt.bar(range(len(dict_sex_count)), dict_sex_count.values(), align="center")
        plt.xticks(range(len(dict_sex_count)), list(map(lambda x: "\n".join(str(x).split(" ")), list(dict_sex_count.keys()))))
        plt.xlabel(col, fontdict={"fontsize":14})
        plt.ylabel("Count", fontdict={"fontsize":14})
        plt.show()
        plt.close()
# %%
