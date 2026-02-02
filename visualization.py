# visualization.py

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os


def plot_population_frequencies(df: pd.DataFrame, output_dir="output"):
    """Bar plot comparing SAS vs EUR vs AFR frequencies for each variant."""

    if df.empty:
        print("No data for plotting.")
        return

    # Melt for seaborn (use your exact column names)
    melted = df.melt(
        id_vars=['CHROM', 'POS', 'END', 'SVTYPE', 'Gene'],
        value_vars=['SAS_AF', 'EUR_AF', 'AFR_AF'],
        var_name='Population',
        value_name='AF'
    )
    melted['AF'] = pd.to_numeric(melted['AF'], errors='coerce').fillna(0)

    plt.figure(figsize=(12, 8))
    sns.barplot(data=melted, x='AF', y='Gene', hue='Population', palette='Set2')
    plt.title("Allele Frequency of SVs in Nutrition Genes by Population")
    plt.xlabel("Allele Frequency")
    plt.ylabel("Gene")
    plt.grid(axis='x', linestyle='--', alpha=0.5)
    plt.legend(title='Population')

    plot_path = os.path.join(output_dir, "sv_population_frequencies.png")
    plt.tight_layout()
    plt.savefig(plot_path, dpi=300)
    plt.close()
    print(f"Population frequency plot saved: {plot_path}")


def plot_sv_size_distribution(df: pd.DataFrame, output_dir="output"):
    """Histogram of SV lengths grouped by gene."""

    if df.empty:
        print("No data for size plot.")
        return

    df['SVLEN_num'] = pd.to_numeric(df['SVLEN'], errors='coerce').abs()
    df = df.dropna(subset=['SVLEN_num'])  # drop invalid

    if len(df) < 3:
        print("Too few variants with valid SVLEN — skipping size plot.")
        return

    plt.figure(figsize=(10, 6))
    sns.histplot(
        data=df,
        x='SVLEN_num',
        hue='Gene',
        multiple='stack',
        bins=30,
        kde=False  # ← Turn off KDE to avoid the "multiple elements" error
    )
    plt.title("Distribution of Structural Variant Sizes by Gene")
    plt.xlabel("Absolute SV Length (bp)")
    plt.ylabel("Count")
    plt.grid(axis='y', linestyle='--', alpha=0.5)

    plot_path = os.path.join(output_dir, "sv_size_distribution.png")
    plt.tight_layout()
    plt.savefig(plot_path, dpi=300)
    plt.close()
    print(f"SV size distribution plot saved: {plot_path}")