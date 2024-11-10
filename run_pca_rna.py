# %%
import glob
import math
import os
import string

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.cm import get_cmap
from pcarna import PCARna

# %%
list_drop_samples = []
pca = PCARna("", "", list_drop_samples=list_drop_samples, list_select_samples=[])
list_study_groups = ["healthy_train", "healthy", "healthy_bgi", "viral_v1", "viral_v2", "viral_v3", "viral_recov", "anxiety", "suicide_attempt", "depression"]
WORKDIR = f"LASSO_INFO/standardized/healthy_illumina/corr_0.3_above_pval_0.05_below/PCA"          
list_meta_columns = ["Sample_Age_Group",
                    "Sample_Phenotype",
                    "Sample_Sex",
                    "Sampling_Year",
                    "Sequencing_Year",
                    "Sequencing_Platform",
                    "Read_Length"]
pcx_num = 1
pcy_num = 2
list_pca_input = glob.glob(f"{WORKDIR}/**/*.tsv", recursive=True)
list_pca_input_filt = list()
for pca_input in list_pca_input:
    x = "_".join(os.path.basename(pca_input).split("_")[3:]).replace(".tsv", "")
    if x in list_study_groups:
        list_pca_input_filt.append(pca_input)
list_pca_exp_input = sorted(list(filter(lambda x: "exp" in x, list_pca_input_filt)))
list_pca_meta_input = sorted(list(filter(lambda x: "meta" in x , list_pca_input_filt)))

list_df_exp_all = list()
list_df_meta_all = list()
for exp, meta in zip(list_pca_exp_input, list_pca_meta_input):
    df_exp = pd.read_csv(exp, sep="\t")
    df_exp_drop_row = df_exp[~df_exp["SUBJID"].isin(list_drop_samples)]
    df_exp_drop_col = df_exp_drop_row.drop(columns=["Sample_Age"])
    list_df_exp_all.append(df_exp_drop_col)

    df_meta = pd.read_csv(meta, sep="\t")
    df_meta = df_meta.rename(columns={"readLength": "Read_Length", "Project_Year": "Sampling_Year", "Sequencing_Date": "Sequencing_Year", "Sequencing_Type_Platform": "Sequencing_Platform", "Sample_Trait": "Sample_Phenotype"})
    df_meta.loc[df_meta["Sample_Phenotype"].str.startswith("Mental"), "Sample_Phenotype"] = "Mental illness"
    df_meta.loc[df_meta["Project-ID"].str.startswith("C19-"), "Sample_Phenotype"] = "Infection"
    df_meta_drop = pca.drop_samples_meta(df_meta)
    df_meta_shortened = pca.shorten_meta_info(df_meta_drop)
    list_df_meta_all.append(df_meta_shortened)

df_exp_all = pd.concat(list_df_exp_all, axis=0)
df_exp_all = df_exp_all.set_index("SUBJID")
df_meta_all = pd.concat(list_df_meta_all, axis=0)
list_samples_for_fit = df_meta_all[df_meta_all["Sample_Phenotype"] == "Healthy"]["Project-ID"].to_list()
list_samples_all = list(df_exp_all.index)

# %%
train_data = "healthy_train"
for exp_train in list_pca_exp_input:
    x = "_".join(os.path.basename(exp_train).split("_")[3:]).replace(".tsv", "")
    if x == train_data:
        df_exp_train = pd.read_csv(exp_train, sep="\t").set_index("SUBJID")
        df_exp_train = df_exp_train.drop(columns=["Sample_Age"])
pca_obj, pca_exp = pca.run_pca(df_exp_train, df_exp_all, std=False)
df_pca_meta = pca.merge_pcadata_and_metadata(pca_exp, df_meta_all)
# %%
cm = 1/2.54
path = "./Arial.ttf"
prop = fm.FontProperties(fname=path)
plt.rcParams['font.family'] = prop.get_name()
plt.rcParams["font.size"] = 5
plot_linewidth = 1
width = 10*cm
height = 16*cm
figsize = (width, height)

fig, axes = pca.draw_pc_biplot(df_pca_meta, 
                               list_meta_columns, 
                               pca_obj, pcx_num, 
                               pcy_num, 
                               dict_legend_ncol = {}, 
                               plot_linewidth = plot_linewidth, 
                               figsize=figsize)
fig.tight_layout(pad=0, h_pad=1.5, w_pad=1.5)    # Adjusted padding to minimize margins while preventing overlap
plt.savefig("Figures/Supplementary_Figure_4.tiff", dpi=600, bbox_inches='tight')
plt.savefig("Figures/Supplementary_Figure_4.png", dpi=600, bbox_inches='tight')
plt.show()
plt.close()

# %%
