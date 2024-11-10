# %%
import os
import pandas as pd

file_qual_encode_table = "/BiO/Store/KOGIC/RNASeq/MetaInfo/Fastq_Base_Quality_Score_Encoding.xlsx"
df_qual_encode_table = pd.read_excel(file_qual_encode_table).astype(str)
dict_qual_encode_table = dict(zip(df_qual_encode_table["Symbol"], df_qual_encode_table["Q-Score"]))
dict_qual_encode_table