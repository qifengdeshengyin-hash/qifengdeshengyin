import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib import cm
from typing import Sequence

class DensityStripePlotter:

    def __init__(
        self,
        infile: str,
        outfile_png: str,
        outfile_rank: str,
        group_order: Sequence[str] = ("pos", "neg", "ref", "mut"),
        bins: int = 5_000,
        cmap_name: str = "Greys",
        vmax: float = 4e-4,
        stripe_height: float = 0.8,
    ):
        self.infile = infile
        self.outfile_png = outfile_png
        self.outfile_rank = outfile_rank
        self.group_order = group_order
        self.bins = bins
        self.cmap_name = cmap_name
        self.vmax = vmax
        self.stripe_height = stripe_height

        self.df = None
        self.N = None
        self.norm = Normalize(vmin=0, vmax=self.vmax)
        self.cmap = cm.ScalarMappable(norm=self.norm, cmap=self.cmap_name).to_rgba

    def load_and_rank(self):
        df = pd.read_csv(self.infile, sep="\t")
        df = df.sort_values("activity", ascending=False).reset_index(drop=True)
        df["rank"] = df.index
        self.df = df
        self.N = len(df)

    def save_rank_table(self):
        self.df.to_csv(self.outfile_rank, sep="\t", index=False)

    def _draw_stripe(self, ax, ranks: np.ndarray, y_pos: float):
        density, bins = np.histogram(ranks, bins=self.bins, range=(0, self.N), density=True)
        segments = [(bins[i], bins[i + 1] - bins[i]) for i in range(len(bins) - 1)]
        colors = [self.cmap(d) for d in density]

        ax.broken_barh(segments, (y_pos - self.stripe_height / 2, self.stripe_height), facecolors=colors, edgecolors="none", )

    def plot(self):
        fig_h = 1.8 * len(self.group_order)
        fig, ax = plt.subplots(figsize=(30, fig_h))

        for idx, grp in enumerate(self.group_order):
            ranks = self.df.loc[self.df["group"] == grp, "rank"].values
            self._draw_stripe(ax, ranks, idx)

        ax.set_yticks(range(len(self.group_order)))
        ax.set_yticklabels(self.group_order, fontsize=18)
        ax.set_xlim(0, self.N)
        ax.set_xlabel("Activity rank", fontsize=18)
        ax.tick_params(axis="x", labelsize=16)
        ax.spines["top"].set_visible(False)
        ax.spines["left"].set_visible(False)

        sm = cm.ScalarMappable(norm=self.norm, cmap=self.cmap_name)
        sm.set_array([])
        cb = plt.colorbar(sm, ax=ax, orientation="vertical", pad=0.01)
        cb.set_label("Density", fontsize=18)
        cb.ax.tick_params(labelsize=14)

        plt.tight_layout()
        plt.savefig(self.outfile_png, dpi=300)
        plt.close()

    def run(self):
        self.load_and_rank()
        self.save_rank_table()
        self.plot()

def main():
    path_in = os.path.expanduser("~/01.SNP-AF/Data/3.data")
    path_out = os.path.expanduser("~/01.SNP-AF/Data/4.plot")

    infile = os.path.join(path_in, "activity_barcode.tsv")
    outfile_png = os.path.join(path_out, "density_stripe_plot_barcode.png")
    outfile_rank = os.path.join(path_out, "rank.tsv")

    plotter = DensityStripePlotter(infile=infile, outfile_png=outfile_png, outfile_rank=outfile_rank, group_order=("pos", "neg", "ref", "mut"), )
    plotter.run()


if __name__ == "__main__":
    main()
