# vcf_parser.py
# Parses 1000 Genomes structural variant VCF (GRCh38)
# Extracts key info: SVTYPE, length, population frequencies (SAS = South Asian)

import os
import pandas as pd
from typing import List, Dict, Optional


def parse_vcf(vcf_path: str,
              gene_regions: Dict = None,  # dict of gene regions from config
              sv_types: List[str] = ["DEL", "DUP"],
              max_variants: Optional[int] = None,
              output_csv: str = "output/sv_results.csv",
              output_txt: str = "output/sv_summary.txt") -> Optional[pd.DataFrame]:
    """
    Parse a VCF file containing structural variants.

    Parameters:
        vcf_path: Path to the unzipped .vcf file
        gene_regions: Optional dict of gene regions (e.g. {"TCF7L2": {"chrom": "10", "start": 112950000, "end": 113300000}})
        sv_types: List of SVTYPE values to keep (e.g. ["DEL", "DUP"])
        max_variants: Limit number of variants processed (for speed during testing)
        output_csv: Where to save structured results as CSV
        output_txt: Where to save human-readable summary

    Returns:
        pandas DataFrame with parsed variants (or None if error)
    """
    if not os.path.exists(vcf_path):
        print(f"Error: VCF file not found at {vcf_path}")
        return None

    print(f"Parsing VCF: {vcf_path}")
    print(f"Looking for SV types: {sv_types}")
    if gene_regions:
        print(f"Filtering for gene regions: {list(gene_regions.keys())}")

    variants: List[Dict] = []
    variant_count = 0
    kept_count = 0

    with open(vcf_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue  # skip empty & header lines

            variant_count += 1
            if max_variants and variant_count > max_variants:
                print(f"Reached max_variants limit ({max_variants})")
                break

            columns = line.split('\t')
            if len(columns) < 9:
                continue

            chrom = columns[0]
            pos = int(columns[1])
            var_id = columns[2]
            ref = columns[3]
            alt = columns[4]
            info = columns[7]

            # Parse INFO field into dict
            info_dict = {}
            for field in info.split(';'):
                if '=' in field:
                    key, value = field.split('=', 1)
                    info_dict[key] = value

            sv_type = info_dict.get('SVTYPE', 'UNKNOWN')
            sv_len = info_dict.get('SVLEN', 'NA')
            end = int(info_dict.get('END', pos))  # fallback to POS if no END
            af_sas = info_dict.get('SAS_AF', 'NA')
            af_eur = info_dict.get('EUR_AF', 'NA')
            af_afr = info_dict.get('AFR_AF', 'NA')

            # Check if variant overlaps any gene region
            matched_gene = None
            if gene_regions:
                for gene, region in gene_regions.items():
                    if chrom == region["chrom"] and pos <= region["end"] and end >= region["start"]:
                        matched_gene = gene
                        break

            # Keep if SV type matches and (optional) gene region matched
            if sv_type in sv_types and (not gene_regions or matched_gene):
                kept_count += 1
                variants.append({
                    'CHROM': chrom,
                    'POS': pos,
                    'END': end,
                    'SVTYPE': sv_type,
                    'SVLEN': sv_len,
                    'SAS_AF': af_sas,
                    'EUR_AF': af_eur,
                    'AFR_AF': af_afr,
                    'Gene': matched_gene or 'None',
                    'INFO': info[:200] + '...' if len(info) > 200 else info
                })

    if not variants:
        print("No matching variants found.")
        with open(output_txt, 'w') as f:
            f.write("No matching structural variants found.\n")
        return None

    # Convert to DataFrame
    df = pd.DataFrame(variants)

    # Create a temporary numeric column for sorting (only if SAS_AF exists)
    if 'SAS_AF' in df.columns:
        df['SAS_AF_num'] = pd.to_numeric(df['SAS_AF'], errors='coerce').fillna(0)
        df = df.sort_values(by='SAS_AF_num', ascending=False)
        df = df.drop(columns='SAS_AF_num')  # clean up temp column

    # Save results
    df.to_csv(output_csv, index=False)
    print(f"Saved {len(df)} variants to {output_csv}")

    # Summary text file
    with open(output_txt, 'w') as f:
        f.write(f"VCF Parsing Summary\n")
        f.write(f"====================\n\n")
        f.write(f"File: {vcf_path}\n")
        f.write(f"Total variants processed: {variant_count}\n")
        f.write(f"Variants kept: {kept_count}\n")
        f.write(f"SV types filtered: {sv_types}\n")
        if gene_regions:
            f.write(f"Gene regions searched: {list(gene_regions.keys())}\n")
        f.write(f"\nFirst 10 rows:\n")
        f.write(df.head(10).to_string(index=False))

    print(f"Summary saved to {output_txt}")

    return df