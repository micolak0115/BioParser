# %%
import os
import pandas as pd

# %%
df_metadata = pd.read_csv("/BiO/Research/GeroPathway/Resources/RawData/Transcriptome/Killifish/GRZ/sra_metadata_run_info_GRZ.txt", sep="\t")
df_tissue = pd.read_csv("/BiO/Research/GeroPathway/Resources/RawData/Transcriptome/Killifish/GRZ/tissue_info_GRZ.txt", sep="\t"
)

list_biosampleid = df_tissue["BioSample_ID"].to_list()
list_biosampleid = list(map(lambda x: f"SAMN{str(x)}", list_biosampleid))
list_tissuename = df_tissue["Tissue"]
dict_sampleid_tissue = dict(zip(list_biosampleid, list_tissuename))

df_metadata["tissue_annot"] = df_metadata["BioSample"].apply(lambda x: dict_sampleid_tissue[x] if x in dict_sampleid_tissue else x)

df_metadata.to_csv("/BiO/Research/GeroPathway/Resources/RawData/Transcriptome/Killifish/GRZ/sra_metadata_run_info_GRZ_tissue_added.txt", sep="\t", index=False)
# %%
