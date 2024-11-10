# %%
import math
import os
import string

import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
import numpy as np
import pandas as pd
from matplotlib.colors import rgb2hex
from matplotlib.patches import Ellipse
from matplotlib.pyplot import get_cmap
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


class PCARna():
    def __init__(self, expmtx, metadata, list_drop_samples=[], list_select_samples=[]):
        self.exptx = expmtx
        self.metadata = metadata
        self.list_drop_samples = list_drop_samples
        self.list_select_samples = list_select_samples
        self.colsampleid = "ID"
        self.sep = "\t"
        self.n_components = 20

    def read_metadata(self):
        df_meta = pd.read_csv(self.metadata, sep=self.sep)

        return df_meta

    def drop_samples_meta(self, df_meta, colsampleid="Project-ID"):
        df_meta_drop = df_meta[~df_meta[colsampleid].isin(self.list_drop_samples)]

        return df_meta_drop

    def shorten_meta_info(self, df_meta):
        df_meta_shortened = df_meta.copy()
        for col in list(df_meta_shortened.columns)[3:]:
            meta_info = df_meta_shortened[col]
            new_meta_info = list()
            for elem in meta_info:
                if len(str(elem)) > 20:
                    new_elem = str(elem)[:21] + ".."
                else:
                    new_elem = str(elem)
                new_meta_info.append(new_elem)
            df_meta_shortened[col] = new_meta_info
        
        return df_meta_shortened

    def read_expmtx(self):
        df_exp = pd.read_csv(self.expmtx, sep=self.sep)
        
        return df_exp

    def drop_samples_exp(self, df_exp):
        df_exp_drop = df_exp.drop(columns=self.list_drop_samples)

        return df_exp_drop

    def select_samples_df(self, df_exp):
        df_exp_indexed = df_exp.set_index(self.colsampleid)
        df_exp_select = df_exp_indexed[self.list_select_samples]
        df_exp_select_reidx = df_exp_select.reset_index(drop=False)

        return df_exp_select_reidx

    def save_expmtx(self, df_save, path_output):
        df_save.to_csv(path_output, sep=self.sep, index=False)

        return None

    def transpose_exp_matrix(self, df_exp, colsampleid=None):
        if colsampleid == None:
            list_geneid = df_exp.index.to_list()
            df_exp_transposed = df_exp.T
            df_exp_transposed.columns = list_geneid

        else:
            list_geneid = df_exp[colsampleid].to_list()
            df_exp_transposed = df_exp.T.iloc[1:, :]
            df_exp_transposed.columns = list_geneid

        return df_exp_transposed

    def run_pca(self, df_train, df_test, std=True):
        train_values = df_train.values
        if std:
            scaler = StandardScaler().fit(df_train)
            train_values = scaler.transform(df_train)
        pca = PCA(n_components=self.n_components).fit(train_values)
        if df_test is None:
            pca_exp = pca.transform(train_values)
        else:
            test_values = df_test.values
            pca_exp = pca.transform(test_values)

        return pca, pca_exp

    def draw_scree_plot(self, pca, outdir, test_name):
        # Draw scree plot
        y = np.cumsum(pca.explained_variance_ratio_)
        x = list(1, range(len(y)), 1)
        plt.plot(x, y)
        plt.xlabel('number of components')
        plt.ylabel('cumulative explained variance')
        plt.show()
        outfigname = os.path.join(outdir, f"pca_scree_plot_{test_name}.png")
        plt.savefig(outfigname, dpi=600)
        plt.close()

    def merge_pcadata_and_metadata(self, pca_exp, df_meta, colsampleid="ID"):
        # merge pcadata and metadata
        df_meta["Sample_Age_Group"] = df_meta["Sample_Age"].apply(lambda x: int(str(x)[0] + "0"))
        df_meta_sample_id_renamed = df_meta.rename(columns={"Project-ID": colsampleid})
        list_samples = df_meta_sample_id_renamed[colsampleid]
        list_PC = list(map(lambda x: f"PC{x}", list(range(1, self.n_components+1, 1))))
        df_pca_exp = pd.DataFrame(pca_exp, columns=list_PC, index=list_samples).reset_index(drop=False).rename(columns={"index":colsampleid})

        df_pca_meta = pd.merge(df_pca_exp, df_meta_sample_id_renamed, how="inner", on=colsampleid)

        return df_pca_meta

    def get_loading_vector(self, df_exp_dropna, pca, outdir, project_info, colsampleid="ID"):
        # Get loading vector
        list_feature_name = df_exp_dropna[colsampleid].to_list()
        list_pc = [np.transpose(pca.components_[i]) for i in range(pca.n_components_)]

        df_load = pd.DataFrame()
        df_load["Gene_Symbol"] = list_feature_name
        for ind, pc in enumerate(list_pc, 1):
            df_load[f"PC{ind}"] = pc

        df_load.to_csv(os.path.join(outdir, f"pca_loading_vector_{project_info}.tsv"))

        return df_load

    def extract_top_pc_features(self, df_exp_transposed, pca, num_top=50):
        # Extract genes involved in making components
        # https://stackoverflow.com/questions/50796024/feature-variable-importance-after-a-pca-analysis
        # name of features
        feature_names = list(df_exp_transposed.columns)
        # number of components
        n_pcs= pca.components_.shape[0]
        # get the index of the most important feature on EACH component
        most_important = [np.abs(pca.components_[i]).argsort()[::-1][:num_top] for i in range(n_pcs)]

        for _, list_idx in enumerate(list(map(list, most_important)), 1):
            list_features_top = list()
            for idx in list_idx:
                feat = feature_names[idx]
                list_features_top.append(feat)
        
        return list_features_top

    def confidence_ellipse(self, x, y, ax, n_std=3.0, facecolor='none', **kwargs):
        """
        Create a plot of the covariance confidence ellipse of *x* and *y*.

        Parameters
        ----------
        x, y : array-like, shape (n, )
            Input data.

        ax : matplotlib.axes.Axes
            The axes object to draw the ellipse into.

        n_std : float
            The number of standard deviations to determine the ellipse's radiuses.

        **kwargs
            Forwarded to `~matplotlib.patches.Ellipse`

        Returns
        -------
        matplotlib.patches.Ellipse
        """
        if x.size != y.size:
            raise ValueError("x and y must be the same size")

        cov = np.cov(x, y)
        pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
        # Using a special case to obtain the eigenvalues of this
        # two-dimensional dataset.
        ell_radius_x = np.sqrt(1 + pearson)
        ell_radius_y = np.sqrt(1 - pearson)
        ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2,
                        facecolor=facecolor, **kwargs)

        # Calculating the standard deviation of x from
        # the squareroot of the variance and multiplying
        # with the given number of standard deviations.
        scale_x = np.sqrt(cov[0, 0]) * n_std
        mean_x = np.mean(x)

        # calculating the standard deviation of y ...
        scale_y = np.sqrt(cov[1, 1]) * n_std
        mean_y = np.mean(y)

        transf = transforms.Affine2D() \
            .rotate_deg(45) \
            .scale(scale_x, scale_y) \
            .translate(mean_x, mean_y)

        ellipse.set_transform(transf + ax.transData)
        return ax.add_patch(ellipse)
    
    def rgb_to_cmyk(self, r, g, b, CMYK_SCALE = 100, RGB_SCALE = 1):
        if (r, g, b) == (0, 0, 0):
            # black
            return 0, 0, 0, CMYK_SCALE

        # rgb [0,255] -> cmy [0,1]
        c = 1 - r / RGB_SCALE
        m = 1 - g / RGB_SCALE
        y = 1 - b / RGB_SCALE

        # extract out k [0, 1]
        min_cmy = min(c, m, y)
        c = (c - min_cmy) / (1 - min_cmy)
        m = (m - min_cmy) / (1 - min_cmy)
        y = (y - min_cmy) / (1 - min_cmy)
        k = min_cmy

        # rescale to the range [0,CMYK_SCALE]
        return c * CMYK_SCALE, m * CMYK_SCALE, y * CMYK_SCALE, k * CMYK_SCALE

    def cmyk_to_rgb(self, c, m, y, k, cmyk_scale = 100, rgb_scale=1, alpha = 1):
        r = rgb_scale * (1.0 - c / float(cmyk_scale)) * (1.0 - k / float(cmyk_scale))
        g = rgb_scale * (1.0 - m / float(cmyk_scale)) * (1.0 - k / float(cmyk_scale))
        b = rgb_scale * (1.0 - y / float(cmyk_scale)) * (1.0 - k / float(cmyk_scale))
        return r, g, b, alpha

    def draw_pc_biplot(self,
                    df_pca_input, 
                    list_meta_columns, 
                    pca_obj, 
                    pcx_num, 
                    pcy_num, 
                    ncol=2, 
                    dict_legend_ncol=dict(),
                    plot_linewidth = 1,
                    figsize=(5, 5)):
        
        pcx = f"PC{pcx_num}"
        pcy = f"PC{pcy_num}"
        nrow = math.ceil(len(list_meta_columns) / ncol)
        fig, axes = plt.subplots(nrows=nrow, ncols=ncol, figsize=figsize)  # Increased figure size
        
        # Flatten axes to handle the case when it's a 2D array
        if nrow > 1 or ncol > 1:
            axes = axes.flatten()
        else:
            axes = np.array([axes]).flatten()
        
        for i, (meta, ax) in enumerate(zip(list_meta_columns, axes)):
            fig_letter = list(string.ascii_uppercase)[i]
            val_groupby_meta_pcx = df_pca_input.groupby(meta)[pcx].apply(np.array).tolist()
            val_groupby_meta_pcy = df_pca_input.groupby(meta)[pcy].apply(np.array).tolist()
            cmap = get_cmap("Reds", len(val_groupby_meta_pcx))
            list_colors = [cmap(ind) for ind in range(len(val_groupby_meta_pcx))]
            dim_level = 10
            list_colors_cmyk = [self.rgb_to_cmyk(rgb[0], rgb[1], rgb[2]) for rgb in list_colors]
            list_colors_cmyk_dim = [(cmyk[0], cmyk[1], cmyk[2], cmyk[3] + dim_level) for cmyk in list_colors_cmyk]
            list_colors_rgb_dim = [self.cmyk_to_rgb(cmyk[0], cmyk[1], cmyk[2], cmyk[3]) for cmyk in list_colors_cmyk_dim]
            list_legends = df_pca_input.groupby(meta)[pcy].apply(np.array).index.tolist()
            list_legends = [int(float(x)) if str(x).endswith(".0") else x for x in list_legends]
            legend_ncol = dict_legend_ncol.get(meta, 1)

            for (x, y, color, legend) in zip(val_groupby_meta_pcx, val_groupby_meta_pcy, list_colors_rgb_dim, list_legends):
                ax.scatter(x, y, c=color, s=plot_linewidth, label=legend, alpha=0.3, zorder=-11)  # Increased marker size
                self.confidence_ellipse(x, y, ax, n_std=3.0, edgecolor=color, facecolor='none', linewidth=plot_linewidth*2.0, alpha=0.5)
                ax.scatter(np.mean(x), np.mean(y), c=color, edgecolor="k", linewidth=plot_linewidth*0.5, marker="*", s=plot_linewidth*70, alpha=1, zorder=999)  # Increased mean marker size
            
            ax.set_xlim(-20, 80)
            ax.set_ylim(-30, 20)
            ax.set_title(meta, 
                        fontsize=plt.rcParams["font.size"]+1,
                        fontweight="bold",
                        pad=2)  # Adjusted padding to avoid overlap with y-tick labels
            ax.set_ylabel(f"principal component {pcy_num}\n({round(pca_obj.explained_variance_ratio_[pcy_num - 1] * 100, 1)}%)", 
                        fontsize=plt.rcParams["font.size"]+1, 
                        labelpad=2)  # Added label padding
            ax.set_xlabel(f"principal component {pcx_num}\n({round(pca_obj.explained_variance_ratio_[pcx_num - 1] * 100, 1)}%)", 
                        fontsize=plt.rcParams["font.size"]+1, 
                        labelpad=2)  # Added label padding
                    
            # Add legend inside the plot without overlapping with data points
            legend = ax.legend(loc="upper right", 
                            fontsize=plt.rcParams["font.size"],
                            ncol=legend_ncol, 
                            markerscale=plot_linewidth, 
                            # borderpad=plot_linewidth-0.5, 
                            handletextpad=1.0, 
                            bbox_to_anchor=(1.0, 1.0), 
                            frameon=True)
            # Set legend box properties
            # legend.set_frame_on(True)  # Ensure the legend box has a frame
            # legend.get_frame().set_facecolor('white')  # Set legend box face color to white
            # legend.get_frame().set_edgecolor('black')  # Set legend box edge color
            # legend.get_frame().set_linewidth(1.5)  # Set legend box edge linewidth
            legend.get_frame().set_alpha(1.0)  # Set legend box alpha to make it opaque
            # Make legend text bold
            plt.setp(legend.get_texts(), weight='bold')
            
            # Make legend dots completely opaque
            for handle in legend.legendHandles:
                handle.set_alpha(1)
            
            # Annotate figure lettering outside the plot, on the far left to avoid overlap
            ax.annotate(fig_letter,
                        xy=(-0.4, 1.04),  # Adjusted further left to avoid overlap
                        xycoords='axes fraction',
                        xytext=(0, 0),
                        textcoords='offset points',
                        size=plt.rcParams["font.size"]+2, 
                        ha='left', 
                        va='center',
                        fontweight="bold",
                        color="black")
            
        for extra_ind in range(len(list_meta_columns), len(axes)):
            axes[extra_ind].axis("off")
        
        return fig, axes

    def get_list_selected_genes(self, infile):
        list_selected_genes = list()
        with open(infile, mode="r") as fr:
            _skiprow = fr.readline()
            for line in fr:
                record = line.rstrip("\n").split("\t")
                gene = record[-1]
                list_selected_genes.append(gene)

        return list_selected_genes

    def save_selected_genes(self, list_selected_genes, outdir):
        fileingenename = os.path.join(outdir, "features_ingenename.txt")
        with open(fileingenename, mode="w") as fw:
            for gene in list_selected_genes:
                fw.write(gene + "\n")

    def get_important_features(self, file_features, colfeat="Selected_Features"):
        list_features = list()
        with open(file_features, mode="r") as fr:
            list_header = fr.readline().rstrip("\n").split("\t")
            idx_feat = list_header.index(colfeat)
            for line in fr:
                record = line.rstrip().split()
                gene = record[idx_feat]
                list_features.append(gene)

        return list_features

    def get_num_original_feature(self, dir_scaled_exp):
        file_original_features = os.path.join(os.path.dirname(dir_scaled_exp), "regression_results.txt")
        list_original_features = self.get_important_features(file_original_features)
        num_original_features = len(list_original_features)

        return num_original_features

    def get_path_file_features(self, dir_scaled_exp):
        list_filename_scaled_exp = os.listdir(dir_scaled_exp)
        list_path_scaled_exp = list(map(lambda x: os.path.join(dir_scaled_exp, x), list_filename_scaled_exp))

        return list_path_scaled_exp

    def get_select_random_features(self, scaled_exp, num_original_features):
        with open(scaled_exp, mode="r") as fr:
            list_header = fr.readline().rstrip("\n").split("\t")
            np.random.seed(1)
            list_rand_features = list(np.random.choice(list_header[1:-1], size=num_original_features, replace=False))
            list_new_header = [list_header[0]] + list_rand_features + [list_header[-1]]

        return list_new_header
# %%
