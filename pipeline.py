import requests
import gzip
import shutil
import os


def setup_output_directory(output_dir="output"):
    """
    Create the output directory if it doesn't exist.
    Returns the full path to the directory.
    """
    os.makedirs(output_dir, exist_ok=True)
    print(f"Output directory ready: {output_dir}")
    return output_dir

def download_sv_vcf(output_dir="output"):
    """Download the 1000 Genomes integrated SV genotypes VCF (GRCh38)."""

    url = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/supporting/GRCh38_positions/ALL.wgs.integrated_sv_map_v2_GRCh38.20130502.svs.genotypes.vcf.gz"

    gz_path = os.path.join(output_dir, "1000g_sv_genotypes_GRCh38.vcf.gz")
    vcf_path = os.path.join(output_dir, "1000g_sv_genotypes_GRCh38.vcf")

    print(f"Downloading SV VCF from {url}...")
    try:
        response = requests.get(url, stream=True, timeout=60)
        response.raise_for_status()

        print(f"HTTP Status: {response.status_code}")
        print(f"Content-Length: {response.headers.get('Content-Length', 'unknown')} bytes")

        with open(gz_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)

        print(f"Downloaded to {gz_path}")

    except requests.exceptions.RequestException as e:
        print(f"Download error: {e}")
        return None

    # Unzip
    print("Unzipping...")
    try:
        with gzip.open(gz_path, 'rb') as f_in:
            with open(vcf_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        print(f"Unzipped to {vcf_path}")
    except Exception as e:
        print(f"Unzip error: {e}")
        return None

    return vcf_path


def explore_vcf(vcf_path, output_file="output/vcf_summary.txt"):
    """Explore VCF file: count lines, extract columns, save to text file."""

    with open(vcf_path, 'r') as f:
        total_lines = 0
        header_lines = 0
        variant_lines = 0
        columns = []

        for line in f:
            total_lines += 1
            if line.startswith('#'):
                header_lines += 1
                if line.startswith('#CHROM'):
                    columns = line.strip().split('\t')
            else:
                variant_lines += 1

        size_mb = os.path.getsize(vcf_path) / (1024 * 1024)  # Size in MB

    # Save to text file
    with open(output_file, 'w') as out:
        out.write("VCF File Summary\n")
        out.write("=" * 20 + "\n")
        out.write(f"File path: {vcf_path}\n")
        out.write(f"Size: {size_mb:.2f} MB\n")
        out.write(f"Total lines: {total_lines}\n")
        out.write(f"Header lines: {header_lines}\n")
        out.write(f"Variant rows: {variant_lines}\n")
        out.write(f"Columns ({len(columns)} total):\n")
        out.write("\n".join(columns) + "\n")

    print(f"Summary saved to {output_file}")

