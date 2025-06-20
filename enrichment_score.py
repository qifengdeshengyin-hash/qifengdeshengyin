import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import random
from statsmodels.stats.multitest import multipletests

class EnrichmentAnalyzer:
    def __init__(self, input_file, output_file, group_labels=('mut', 'ref'), n_permutations=10_000):
        self.input_file = input_file
        self.output_file = output_file
        self.group_labels = group_labels
        self.n_permutations = n_permutations
        self.df = None
        self.N = None
        self.group_flags = {}
        self.group_scores = {}
        self.group_means = {}
        self.group_pvals = {}
        self.group_pvals_bonf = {}

    def load_and_prepare(self):
        df = pd.read_csv(self.input_file, sep='\t')
        df = df.sort_values(by='activity', ascending=False).reset_index(drop=True)
        self.df = df
        self.N = len(df)
        for group in self.group_labels:
            self.df[f'is_{group}'] = self.df['group'] == group
            self.group_flags[group] = self.df[f'is_{group}']

    def compute_enrichment_score(self, flags: pd.Series, n_total: int) -> np.ndarray:
        scores, hit = [], 0
        for rank in range(1, self.N + 1):
            if flags[rank - 1]:
                hit += 1
            actual = hit / n_total if n_total else 0
            expected = rank / self.N
            scores.append(actual - expected)
        return np.asarray(scores)

    def permutation_pvalue(self, observed_mean: float, n_true: int) -> float:
        hits = 0
        for _ in range(self.n_permutations):
            rand_idx = random.sample(range(self.N), n_true)
            flags = np.zeros(self.N, dtype=bool)
            flags[rand_idx] = True
            rand_mean = self.compute_enrichment_score(flags, n_true).mean()
            if observed_mean > 0 and rand_mean >= observed_mean:
                hits += 1
            elif observed_mean < 0 and rand_mean <= observed_mean:
                hits += 1
        return hits / self.n_permutations

    def run_analysis(self):
        for group in self.group_labels:
            flags = self.group_flags[group]
            n_group = flags.sum()
            scores = self.compute_enrichment_score(flags, n_group)
            mean_score = scores.mean()
            pval = self.permutation_pvalue(mean_score, n_group)
            self.group_scores[group] = scores
            self.group_means[group] = mean_score
            self.group_pvals[group] = pval

        pvals = list(self.group_pvals.values())
        pvals_bonf = multipletests(pvals, method='bonferroni')[1]
        self.group_pvals_bonf = dict(zip(self.group_labels, pvals_bonf))

    def plot(self):
        plt.figure(figsize=(8, 4))
        for group in self.group_labels:
            scores = self.group_scores[group]
            p_bonf = self.group_pvals_bonf[group]
            color = 'brown' if group == 'mut' else 'blue'
            label = f'{group} (Bonf. p={p_bonf:.3g})'
            plt.plot(range(1, self.N + 1), scores, label=label, color=color)
        plt.axhline(0, color='gray', linestyle='--', linewidth=0.5)
        ax = plt.gca()
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        plt.xlabel('activity rank')
        plt.ylabel('Enrichment score')
        plt.legend(frameon=False)
        plt.tight_layout()
        plt.savefig(self.output_file, dpi=300)
        plt.close()

def main():
    path_in  = os.path.expanduser('~/01.SNP-AF/Data/3.data')
    path_out = os.path.expanduser('~/01.SNP-AF/Data/4.plot')

    file_in  = os.path.join(path_in,  'activity_barcode.tsv')
    # file_out = os.path.join(path_out, 'enrichment_curve_barcode_mut_ref.png')#barcode #region
    file_out = os.path.join(path_out, 'enrichment_curve_barcode_mut_ref_inall.png')

    analyzer = EnrichmentAnalyzer(input_file=file_in, output_file=file_out, group_labels=('mut', 'ref'),  n_permutations=10000)

    analyzer.load_and_prepare()
    analyzer.run_analysis()
    analyzer.plot()

if __name__ == '__main__':
    main()
