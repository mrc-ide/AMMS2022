import os
import click
import uuid
import subprocess
import pandas as pd
from dataclasses import dataclass

# Notes
# - This script subsamples SNPs and samples from MalariaGEN Pf3k data
# - You need to have downloaded the data locally (https://www.malariagen.net/data/pf3K-5)
# - To run, see `python sample_from_pf3k.py --help`
# - Omissions include:
#  - Not filtering by linkage (`bcftools +prune`)
#  - Not filterinng by missingess (`bcftools fiter`)
# JHendry, 2022/07/07


# --------------------------------------------------------------------------------
# SNP SET DEFINITIONS
#
# --------------------------------------------------------------------------------


@dataclass
class SNPSet:
    country: str
    min_af: float
    max_af: float
    min_dp: int = 10
    max_dp: int = 150
    n_snps: int = 96


SNP_SET_COLLECTION = {
    "setA": SNPSet("DRCongo", min_af=0.1, max_af=0.9),
    "setB": SNPSet("DRCongo", min_af=0.3, max_af=0.6),
    "setC": SNPSet("DRCongo", min_af=0.01, max_af=0.05),
    "setD": SNPSet("Vietnam", min_af=0.1, max_af=0.9),
}

# --------------------------------------------------------------------------------
# GLOBAL CONSTANTS
#
# --------------------------------------------------------------------------------


DATA_DIR = "/Users/jasongms/Documents/phd/projects/pf-RMCL/data/pf3k"
METADATA_PATH = f"{DATA_DIR}/metadata/pf3k_release_5_metadata_20190910_cleaned.csv"
VCF_DIR = f"{DATA_DIR}/vcfs/"
VCF_FN_TEMPLATE = "SNP_INDEL_Pf3D7_{chrom:02d}_v3.high_quality_biallelic_snps.vcf.gz"
CHROMS = list(range(1, 15))
OUTPUT_DIR = "../data"


# --------------------------------------------------------------------------------
# SAMPLE SELECTION
#
# --------------------------------------------------------------------------------


def get_country_sample_names(metadata_path, country):
    """
    Get sample IDs for a target country

    """
    # Load
    df = pd.read_csv(metadata_path)
    assert "country" in df.columns
    assert "sample" in df.columns

    # Query
    df.query("country == @country", inplace=True)
    assert df.shape[0] > 0, f"No samples found for country {country}."

    return df["sample"].values.tolist()


# --------------------------------------------------------------------------------
# BCFTOOLS INTERFACE
#
# --------------------------------------------------------------------------------


def bcftools(subcommand, args, input_vcf, run=False):
    """
    Construct `bcftools` commands from python

    """
    assert isinstance(args, list)
    cmd = "bcftools"
    cmd += f" {subcommand}"
    cmd += f" {' '.join([str(a) for a in args])}"
    cmd += f" {input_vcf}"

    if not run:
        return cmd

    subprocess.run(cmd, shell=True, check=True)


# --------------------------------------------------------------------------------
# DOWNSAMPLING
#
# --------------------------------------------------------------------------------


def randomly_downsample_snps_from_vcf(input_vcf, output_vcf, n_snps):
    """
    Randomly downsample to `n_snps` from an `input_vcf`
    (surprisingly can't find bcftools subcommand to do this)

    """
    temp_vcf = f"temp.{str(uuid.uuid4())}.vcf"

    # First, store the header
    cmd = f"bcftools view -h {input_vcf} -o {temp_vcf}"
    subprocess.run(cmd, shell=True, check=True)

    # Randomly downsample by shuffling
    cmd = f"bcftools view -H {input_vcf}"  # view without header
    cmd += f" | shuf -n {n_snps}"  # shuffle and return `n_snps`
    cmd += f" >> {temp_vcf}"  # append to the `output_vcf`
    subprocess.run(cmd, shell=True, check=True)

    # Sort the result
    cmd = f"bcftools sort {temp_vcf} -o {output_vcf}"
    subprocess.run(cmd, shell=True, check=True)
    os.remove(temp_vcf)


# --------------------------------------------------------------------------------
# MAIN SCRIPT
#
# --------------------------------------------------------------------------------


@click.command()
@click.option(
    "-s",
    "--snp_set_name",
    type=click.Choice(SNP_SET_COLLECTION),
    required=True,
    help="Choose SNP set by name. Different sets have differet filters. See script for details.",
)
def main(snp_set_name):
    """
    Grab a reduced set of samples and SNPs from MalariaGEN Pf3k data

    Samples are selected by country
    SNPs are first PLAF and depth filtered, then randomly downsampled.

    """

    # EXTRACT SAMPLES
    snp_set = SNP_SET_COLLECTION[snp_set_name]
    samples = get_country_sample_names(METADATA_PATH, snp_set.country)

    # FILTER FOR EACH CHROMOSOME
    output_vcfs = []
    for chrom in CHROMS:
        input_vcf = f"{VCF_DIR}/{VCF_FN_TEMPLATE.format(chrom=chrom)}"
        output_vcf = os.path.basename(input_vcf).replace(
            ".vcf", f".{snp_set.country}.{snp_set_name}.filtered.vcf"
        )
        output_vcfs.append(output_vcf)

        cmd_samples = bcftools(
            subcommand="view", args=["-s", ",".join(samples)], input_vcf=input_vcf
        )

        cmd_plaf = bcftools(
            subcommand="view",
            args=["--min-af", snp_set.min_af, "--max-af", snp_set.max_af],
            input_vcf="-",
        )

        cmd_depth = bcftools(
            subcommand="view",
            args=[
                "-i",
                f"'MEAN(FORMAT/DP) > {snp_set.min_dp} & MEAN(FORMAT/DP) < {snp_set.max_dp}'",
            ],
            input_vcf="-",
        )

        cmd_clean = bcftools(
            subcommand="annotate",
            args=["-x", "INFO,^FORMAT/GT,FORMAT/AD"],
            input_vcf="-",
        )

        cmd = f"{cmd_samples} | {cmd_plaf} | {cmd_depth} | {cmd_clean} -o {output_vcf}"
        subprocess.run(cmd, shell=True, check=True)

    # CONCATENATE ACROSS CHROMOSOMES
    output_all_vcf = f"SNPs.{snp_set.country}.{snp_set_name}.all_snps.vcf"
    output_sample_vcf = output_all_vcf.replace("all_snps", f"random{snp_set.n_snps}")
    bcftools(
        subcommand="concat",
        args=["-o", output_all_vcf],
        input_vcf=" ".join(output_vcfs),
        run=True,
    )

    # Clean-up post-concateation
    for vcf in output_vcfs:
        os.remove(vcf)

    # DOWNSAMPLE
    randomly_downsample_snps_from_vcf(
        input_vcf=output_all_vcf,
        output_vcf=f"{OUTPUT_DIR}/{output_sample_vcf}",
        n_snps=snp_set.n_snps,
    )


if __name__ == "__main__":
    main()
