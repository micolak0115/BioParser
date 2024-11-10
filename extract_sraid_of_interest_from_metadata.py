# %%
import os
import pandas as pd
import numpy as np
from collections import Counter

# %%
path_metadata = "/BiO/Research/GeroPathway/Resources/RawData/Transcriptome/Killifish/GRZ/sra_metadata_run_info_GRZ_tissue_added.txt"
df_metadata = pd.read_csv(path_metadata, sep="\t")
cond_rna = df_metadata["LibraryStrategy"] == "RNA-Seq"
cond_paired = df_metadata["LibraryLayout"] == "PAIRED"
cond_unknown_tissue = np.logical_or(df_metadata["tissue_annot"].str.contains("SAM"), df_metadata["tissue_annot"].str.contains(","), df_metadata["tissue_annot"].str.contains("whole"))

df_metadata_rna = df_metadata[cond_rna]
df_metadata_rna_paired = df_metadata_rna[cond_paired]
df_metadata_rna_paired_tissue = df_metadata_rna_paired[~cond_unknown_tissue]
# dict_tissue_cnt = Counter(df_metadata_rna_paired_tissue["tissue_annot"])
list_run_id = df_metadata_rna_paired_tissue["Run"].to_list()
list_run_id
# %%
