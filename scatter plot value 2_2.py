import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
    
class FPMAnalyzer:
    def __init__(self, group=None, base_dir="~/01.SNP-AF/Data"):
        self.group = group
        self.base_dir = os.path.expanduser(base_dir)
        self.path1 = os.path.join(self.base_dir, "3.data")
        self.path2 = os.path.join(self.base_dir, "4.plot/normalization")
        # self.DNA_file = os.path.join(self.path1, f"DNA_FPM_matrix{f'_{group}' if group else ''}.tsv")
        self.RNA_file = os.path.join(self.path1, f"RNA_FPM_matrix{f'_{group}' if group else ''}.tsv")
        self.DNA_file = os.path.join(self.path1, f"DNA_FPM_power0.36{f'_{group}' if group else ''}.tsv")
        self.metrics_df = None
        self.region_names = None
        self.annotated_file = os.path.join(self.path2, f"RNA_FPM_matrix_with_group{f'_{self.group}' if self.group else ''}.tsv")

    def load_and_filter_data(self, threshold=0):
        RNA_df = pd.read_csv(self.RNA_file, sep="\t", index_col=0)
        DNA_df = pd.read_csv(self.DNA_file, sep="\t", index_col=0)
        RNA_filtered = RNA_df[(RNA_df > threshold).all(axis=1)]
        DNA_filtered = DNA_df[(DNA_df > threshold).all(axis=1)]
        common_regions = RNA_filtered.index.intersection(DNA_filtered.index)
        self.RNA_df = RNA_filtered.loc[common_regions].copy()
        self.DNA_df = DNA_filtered.loc[common_regions].copy()
        self.region_names = common_regions.tolist()

    def compute_metrics(self):
        RNA_mean = self.RNA_df.mean(axis=1)
        DNA_mean = self.DNA_df.mean(axis=1)
        self.metrics_df = pd.DataFrame({
            "RNA mean": RNA_mean,
            "DNA mean": DNA_mean,
            "RNA / DNA": RNA_mean / DNA_mean,
            # "(RNA+1) / (DNA+1)": (RNA_mean + 1) / (DNA_mean + 1),
            # "log2(RNA / DNA)": np.log2(RNA_mean / DNA_mean),
            # "log2(DNA+1)": np.log2((DNA_mean + 1)),
            # "log2((RNA+1) / (DNA+1))": np.log2((RNA_mean + 1) / (DNA_mean + 1)),
            # "log2(RNA mean)": np.log2(RNA_mean),
            # "log2(DNA mean)": np.log2(DNA_mean),
            # "(RNA+1)": (RNA_mean + 1),
            # "(DNA+1)**0.3": ((DNA_mean + 1)**0.3),
            # "sqrt(DNA+1)": (np.sqrt(DNA_mean + 1)),
            # "(DNA+1)**0.7": ((DNA_mean + 1)**0.7),
            # "ln(DNA+1)": np.log(DNA_mean + 1),
            # "log2(RNA**0.9 / DNA**0.3)": np.log2(RNA_mean**0.9 / DNA_mean**0.3),
            # "log10(DNA+1)": np.log10(DNA_mean + 1),
            # "(RNA+1)**0.3 / (DNA+1)**0.3": ((RNA_mean + 1)**0.3 / (DNA_mean + 1)**0.3),
            # "sqrt(RNA+1) / sqrt(DNA+1)": (np.sqrt(RNA_mean + 1) / np.sqrt(DNA_mean + 1)),
            # "(RNA+1)**0.7 / (DNA+1)**0.7": ((RNA_mean + 1)**0.7 / (DNA_mean + 1)**0.7),
            # "(RNA+1)**0.7 / (DNA+1)**0.3": ((RNA_mean + 1)**0.7 / (DNA_mean + 1)**0.3),
            # "(RNA+1)**0.3 / (DNA+1)**0.7": ((RNA_mean + 1)**0.3 / (DNA_mean + 1)**0.7),
            # "log2(RNA+1) / log2(DNA+1)": np.log2(RNA_mean + 1) / np.log2(DNA_mean + 1),
            # "ln(RNA+1) / ln(DNA+1)": np.log(RNA_mean + 1) / np.log(DNA_mean + 1),
            # "log10(RNA+1) / log10(DNA+1)": np.log10(RNA_mean + 1) / np.log10(DNA_mean + 1),
            # "(RNA+1) / log2(DNA+1)": (RNA_mean + 1) / np.log2(DNA_mean + 1),
        }, index=self.RNA_df.index)

    def annotate_groups(self):
        self.metrics_df['region_barcode'] = self.metrics_df.index
        self.metrics_df['group'] = self.metrics_df['region_barcode'].apply(
            lambda x: x[-4:-1] if x.startswith('rs') else x[:3] if x.startswith(('pos', 'neg')) else 'unknown'
        )
        self.metrics_df.to_csv(self.annotated_file, sep="\t")

    def safe_filename(self, s):
        return (
            s.replace(" ", "").replace("/", "_over_").replace("+", "p")
             .replace("(", "").replace(")", "").replace("sqrt", "sqrt")
             .replace("log2", "log2").replace("*", "x").replace("=", "")
             .replace(",", "")
        )

    def scatter_plot_subplots(self, df, y_key, x_key, xlim, ylim, group_colors, output_prefix):
        fig, axes = plt.subplots(2, 2, figsize=(30, 20))
        highlight_groups = list(group_colors.keys())

        for ax, grp in zip(axes.flatten(), highlight_groups):
            for group in df['group'].unique():
                sub_df = df[df['group'] == group]
                color = group_colors[group] if group == grp else '#404040'
                alpha_val = 0.9 if group == grp else 0.3
                sns.scatterplot(
                    data=sub_df,
                    x=x_key, y=y_key,
                    alpha=alpha_val, s=10, edgecolor='none', color=color, ax=ax, label=group if group == grp else None
                )

            ax.set_title(f"powerDNA0.36_RNA: {grp}", fontsize=18)
            ax.set_xlabel(x_key, fontsize=18)
            ax.set_ylabel(y_key, fontsize=18)
            ax.grid(True, linestyle='--', alpha=0.9)
            # if xlim and ylim:
            #     pad_x = 0.8 * (xlim[1] - xlim[0])
            #     pad_y = 0.8 * (ylim[1] - ylim[0])
            #     ax.set_xlim(xlim[0] - pad_x, xlim[1] + pad_x)
            #     ax.set_ylim(ylim[0] - pad_y, ylim[1] + pad_y)

        plt.tight_layout()
        file_name = f"{output_prefix}{self.safe_filename(y_key)}_vs_{self.safe_filename(x_key)}_4in1.pdf"
        output_file = os.path.join(self.path2, file_name)
        fig.savefig(output_file, dpi=300)
        plt.close(fig)

    def run(self):
        self.load_and_filter_data()
        self.compute_metrics()
        self.annotate_groups()

        df = pd.read_csv(self.annotated_file, sep="\t", index_col=0)
        df['region_barcode'] = df.index

        group_colors = {'pos': 'orangered', 'neg': 'cyan', 'mut': 'lightgreen', 'ref': 'plum'}

        plot_pairs = [
            ("RNA mean", "DNA mean", (0, 700), (0, 350)),
            ("RNA / DNA", "DNA mean", (0, 700), (0, 5.5)),
            # ("log2(RNA**0.9 / DNA**0.3)", "log2(DNA mean)", (3.5, 9.5), (5.0, 9.0)),
            # ("(RNA+1) / (DNA+1)", "DNA mean", (0, 650), (0, 5.5)),
            # ("log2(RNA mean)", "log2(DNA mean)", (3.5, 9.5), (5.0, 9.0)),
            # ("log2(RNA / DNA)", "log2(DNA mean)", (3.5, 10), (-2.5, 2.5)),
            # ("log2((RNA+1) / (DNA+1))", "log2(DNA+1)", (3.5, 10), (-2.5, 2.5)),
            # ("(RNA+1)**0.3 / (DNA+1)**0.3", "(DNA+1)**0.3", (3.5, 9.5), (1.5, 4)),
            # ("(RNA+1)**0.3 / (DNA+1)**0.3", "log2(DNA+1)", (3.5, 9.5), (1.5, 4)),
            # ("sqrt(RNA+1) / sqrt(DNA+1)", "sqrt(DNA+1)", (3.5, 10), (-2.5, 2.5)),
            # ("sqrt(RNA+1) / sqrt(DNA+1)", "log2(DNA+1)", (3.5, 10), (-2.5, 2.5)),
            # ("(RNA+1)**0.7 / (DNA+1)**0.7", "(DNA+1)**0.7", (3.5, 9.5), (3.5, 7)),
            # ("(RNA+1)**0.7 / (DNA+1)**0.7", "log2(DNA+1)", (3.5, 9.5), (3.5, 7)),
            # ("(RNA+1)**0.7 / (DNA+1)**0.3", "(DNA+1)**0.3", (3.5, 10), (-2.5, 2.5)),
            # ("(RNA+1)**0.7 / (DNA+1)**0.3", "log2(DNA+1)", (3.5, 10), (-2.5, 2.5)),
            # ("(RNA+1)**0.3 / (DNA+1)**0.7", "(DNA+1)**0.7", (3.5, 10), (-2.5, 2.5)),
            # ("log2(RNA+1) / log2(DNA+1)", "log2(DNA+1)", (3.5, 10), (-2.5, 2.5)),
            # ("(RNA+1) / log2(DNA+1)", "log2(DNA+1)", (3.5, 10), (-2.5, 2.5)),
            # ("ln(RNA+1) / ln(DNA+1)", "ln(DNA+1)", (3.5, 10), (-2.5, 2.5)),
            # ("ln(RNA+1) / ln(DNA+1)", "log2(DNA+1)", (3.5, 10), (-2.5, 2.5)),
            # ("log10(RNA+1) / log10(DNA+1)", "log10(DNA+1)", (3.5, 10), (-2.5, 2.5)),
            # ("log10(RNA+1) / log10(DNA+1)", "log2(DNA+1)", (3.5, 10), (-2.5, 2.5)),
        ]

        for y_key, x_key, xlim, ylim in plot_pairs:
            self.scatter_plot_subplots(df, y_key, x_key, xlim, ylim, group_colors, output_prefix="powerDNA0.36_RNA")

if __name__ == "__main__":
    analyzer = FPMAnalyzer()
    analyzer.run()
