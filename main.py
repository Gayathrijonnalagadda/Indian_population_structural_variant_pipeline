from pipeline import download_sv_vcf, setup_output_directory, explore_vcf
from vcf_parser import parse_vcf
from config import GENE_REGIONS
from visualization import plot_population_frequencies, plot_sv_size_distribution


def main():
    print("Starting Indian Population Genetic Variation Pipeline\n")

    output_dir = setup_output_directory("output")

    vcf_path = download_sv_vcf(output_dir)

    if vcf_path:
        print(f"\nVCF successfully downloaded and unzipped: {vcf_path}")

        print("\nExploring VCF file...")
        explore_vcf(vcf_path)

        print("\nParsing VCF for gene regions and SV types...")
        df = parse_vcf(
            vcf_path,
            gene_regions=GENE_REGIONS,
            sv_types=["DEL", "DUP"],
            max_variants=None  # Full file
        )

        if df is not None and not df.empty:
            print(f"Found {len(df)} matching variants!")
            print(df.head(10))  # Preview in console

            # Plots — safe inside this if
            plot_population_frequencies(df)
            plot_sv_size_distribution(df)
        else:
            print("No matching variants found or parsing failed.")
    else:
        print("Download failed — check URL or internet connection.")


if __name__ == "__main__":
    main()