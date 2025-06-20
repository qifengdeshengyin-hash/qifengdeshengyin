import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns

path1 = os.path.expanduser('~/01.SNP-AF/Data/3.data')
path2 = os.path.expanduser('~/01.SNP-AF/Data/4.plot')

class CoveragePlotter:
    def __init__(self, count_file, distribution_plot):
        self.input_path = os.path.join(path1, count_file)
        self.output_path = os.path.join(path2, distribution_plot)

    def plot_distribution(self):
        df = pd.read_csv(self.input_path, sep='\t')
        sns.set_style("whitegrid", {'axes.edgecolor': '0.8'})
        plt.rcParams.update({'font.size':12, 'axes.labelsize':14, 'axes.titlesize':16, 'xtick.labelsize':12,'ytick.labelsize':12})

        fig, ax = plt.subplots(figsize=(8, 6))
        counts, bins, _ = ax.hist(df['FPM'], bins=150, color='#4C72B0', edgecolor='black', alpha=0.85, linewidth=1)

        ax.axvline(x=20, color='red', linestyle='-', linewidth=1)
        ymax = counts.max()
        ax.annotate('Coverage > 20\n8082 / 8960 regions', xy=(30, ymax*0.92), xytext=(200, ymax*0.92),
            arrowprops=dict(arrowstyle='->', color='red', linewidth=1 ),
            fontsize=12, color='black', ha='left', va='center')

        ax.set_xlabel('Read coverage (FPM)', labelpad=10)
        ax.set_ylabel('Number of regions', labelpad=10)
        sns.despine(ax=ax) 
        ax.grid(None)
        ax.spines['bottom'].set_color('black')
        ax.spines['left'].set_color('black')
        ax.set_ylim(-0.04*ymax, ymax*1.05)
        ax.margins(x=0.04) 
        ax.tick_params(axis='both', which='major', direction='out', length=3, width=1, colors='black', bottom=True, top=False, left=True, right=False)
        fig.tight_layout(pad=3)
        plt.savefig(self.output_path, dpi=600, bbox_inches='tight')
        plt.close()


def main():
    plotter = CoveragePlotter(count_file='DNA_FPM_mean.tsv', distribution_plot='DNA_FPM_mean.png')
    # plotter = CoveragePlotter(count_file='RNA_FPM_mean.tsv', distribution_plot='RNA_FPM_mean.png')
    # plotter = CoveragePlotter(count_file='DNA_FPM_power0.36.tsv', distribution_plot='DNA_FPM_power0.36.png')
    plotter.plot_distribution()

if __name__ == "__main__":
    main()

