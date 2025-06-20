import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from itertools import combinations

class CorrelationPlotter:
    def __init__(self, input_dir, output_dir):
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.file_mapping = {"DNA": "DNA_FPM_matrix.tsv", "RNA": "RNA_FPM_matrix.tsv", "Activity": "log2_activity_matrix.tsv"}


    def load_data(self, filename):
        file_path = os.path.join(self.input_dir, filename)
        df = pd.read_csv(file_path, sep='\t')
        df.columns = ['region_barcode', 'rep1', 'rep2', 'rep3']
        return df[['rep1', 'rep2', 'rep3']]

    def compute_corr(self, df):
        return df.corr()

    def plot_correlation_matrix(self, df, title, outfile, label):
        corr = df.corr()
        pairs = list(combinations(df.columns, 2))

        fig, axes = plt.subplots(3, 3, figsize=(6, 6))
        color_map = {
            "DNA": "#0D46A1", "RNA": "#096a09", "Activity": "#A76009"}

        for i in range(3):
            for j in range(3):
                ax = axes[i, j]
                if i == j:
                    ax.text(0.5, 0.5, f"A{i+1}", fontsize=20, fontweight='bold', color=color_map[label], ha='center', va='center')
                elif i < j:
                    ax.scatter(df.iloc[:, j], df.iloc[:, i], s=2, color='black')
                    ax.margins(x=0.05, y=0.05)
                    ax.tick_params(axis='both', which='both',direction='out',
                        length=3, width=0.8,labelsize=8,top=True, right=True,)
                    ax.xaxis.set_ticks_position('both')
                    ax.yaxis.set_ticks_position('both')
                    ax.tick_params(labeltop=True, labelright=True)
                else:
                    r = corr.iloc[i, j]
                    ax.text(0.5, 0.5, f"{r:.2f}", fontsize=16, ha='center', va='center', color=color_map[label])

                ax.set_xticks([])
                ax.set_yticks([])
                for spine in ax.spines.values():
                    spine.set_color('gray')
                    spine.set_linewidth(0.5)

        plt.tight_layout()
        plt.suptitle(title, fontsize=14, y=1.02)  
        fig.savefig(outfile, dpi=300, bbox_inches='tight')
        plt.close()

    def run(self):
        for label, filename in self.file_mapping.items():
            df = self.load_data(filename) 
            outfile = os.path.join(self.output_dir, f"{label}_replicate_correlation.png")
            self.plot_correlation_matrix(df,
                title=f"{label} replicate correlation", outfile=outfile, label=label)

def main():
    path1 = os.path.expanduser('~/01.SNP-AF/Data/3.data')
    path2 = os.path.expanduser('~/01.SNP-AF/Data/4.plot')
    plotter = CorrelationPlotter(path1, path2)
    plotter.run()

if __name__ == "__main__":
    main()
