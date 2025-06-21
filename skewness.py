import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import skew as scipy_skew

class SkewnessNormalizer:
    def __init__(self, input_path, output_path, plot_path):
        self.input_path = input_path
        self.output_path = output_path
        self.plot_path = plot_path
        os.makedirs(self.output_path, exist_ok=True)
        os.makedirs(self.plot_path, exist_ok=True)
        self.data = None

    def load_data(self, filename):
        file_path = os.path.join(self.input_path, filename)
        df = pd.read_csv(file_path, sep='\t')
        df.columns = ['region_barcode', 'FPM']
        self.data = df
        return df

    def calculate_skew(self, values):
        return scipy_skew(values, bias=True)

    def apply_power_transform(self, power):
        transformed = np.power(self.data['FPM'].values + 1e-6, power)
        return transformed
    def find_power_for_target_skewness(self, target=1.0, power_range=(0.1, 1.0), step=0.01, tol=0.01, log_file='power_skew_log.tsv'):
        best_p = None
        best_skew = None
        best_data = None
        log_records = []

        powers = np.arange(power_range[1], power_range[0] - step, -step)  # 倒序查找

        for p in powers:
            transformed = self.apply_power_transform(p)
            skew = self.calculate_skew(transformed)

            log_records.append({'power': round(p, 4), 'skew': round(skew, 6)})

            if best_skew is None or abs(skew - target) < abs(best_skew - target):
                best_p = p
                best_skew = skew
                best_data = transformed

        log_df = pd.DataFrame(log_records)
        log_df.to_csv(os.path.join(self.output_path, log_file), sep='\t', index=False)

        return best_p, best_data, best_skew

    def save_transformed_data(self, transformed_values, power):
        df_out = pd.DataFrame({
            'region_barcode': self.data['region_barcode'],
            'FPM': transformed_values
        })
        filename = f'DNA_FPM_power{round(power, 2)}.tsv'
        df_out.to_csv(os.path.join(self.output_path, filename), sep='\t', index=False)

    def plot_distribution(self, values, power):
        plt.figure()
        plt.hist(values, bins=100, edgecolor='black')
        plt.xlabel('FPM')
        plt.ylabel('Counts')
        plt.title(f'Power Transformed Distribution (p={round(power, 2)})')
        plt.tight_layout()
        filename = f'DNA_FPM_power{round(power, 2)}_distribution.png'
        plt.savefig(os.path.join(self.plot_path, filename),dpi=300)
        plt.close()

def main():
    path1 = os.path.expanduser('~/01.SNP-AF/Data/3.data')
    path2 = os.path.expanduser('~/01.SNP-AF/Data/4.plot')

    normalizer = SkewnessNormalizer(input_path=path1, output_path=path1, plot_path=path2)
    df = normalizer.load_data('DNA_FPM_mean.tsv')
    print(f"Loaded {len(df)} rows.")

    power, transformed, final_skew = normalizer.find_power_for_target_skewness(target=1.0, step=0.01, tol=0.01)
    if power is not None:
        print(f"Found power: {power:.2f} with skewness: {final_skew:.4f}")
        print(f"Found power: {power:.2f}")
        normalizer.save_transformed_data(transformed, power)
        normalizer.plot_distribution(transformed, power)
        print(f"Saved transformed data and plot for power={power:.2f}")
    else:
        print("No power found that achieves skewness = 1. Try increasing tolerance or power range.")

if __name__ == '__main__':
    main()
