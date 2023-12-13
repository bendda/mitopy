from .executable import Executable
from .utils import (
    get_file_basename,
    get_file_directory,
    create_output_path,
    check_files_exist,
)
import logging
import os
import pandas as pd
import sys
import glob
from .constants import MT_REFS, MT_BLACKLIST


def _get_contamination(vcf: str, out_dir: str, haplocheck_exec: Executable) -> float:
    """Estimate sample contamination using haplocheck."""

    out_fn = f"{out_dir}/tmp_haplo"
    params = {"--raw": True, "--out": out_fn}

    logging.info("Estimating sample contamination level using Haplocheck...")
    haplocheck_exec.run(vcf, **params)

    contamination_estimate = pd.read_csv(f"{out_fn}.raw.txt", sep="\t")[
        "Contamination Level"
    ].values[0]

    for f in glob.glob(f"{out_fn}*"):
        os.remove(f)
    return 0.0 if contamination_estimate == "ND" else contamination_estimate


def _filter_by_params(
    vcf: str,
    mt_ref: str,
    stats: str,
    vaf_treshold: float,
    max_alt_allele_count: int,
    f_score_beta: float,
    contamination: float,
    out_fn: str,
    gatk_exec: Executable,
) -> None:
    """Filter variants by parameters using gatk FilterMutectCalls."""

    params = {
        "-V": vcf,
        "-R": mt_ref,
        "-O": out_fn,
        "--stats": stats,
        "--max-alt-allele-count": max_alt_allele_count,
        "--mitochondria-mode": True,
        "--min-allele-fraction": vaf_treshold,
        "--f-score-beta": f_score_beta,
        "--contamination-estimate": contamination,
    }

    gatk_exec.run(subcommand="FilterMutectCalls", **params)


def _filter_blacklist(
    vcf: str, out_fn: str, blacklisted_sites: str, gatk_exec: Executable
) -> None:
    """Filter out blacklisted sites using gatk VariantFiltration."""

    params = {
        "-V": vcf,
        "-O": out_fn,
        "--apply-allele-specific-filters": True,
        "--mask": blacklisted_sites,
        "--mask-name": "blacklisted_site",
    }

    gatk_exec.run(subcommand="VariantFiltration", **params)


def _filter_numts(
    vcf: str, mt_ref: str, out_fn: str, autosomal_coverage: str, gatk_exec: Executable
) -> None:
    """Filter out NuMTs using gakt NuMTFilterTool."""

    params = {
        "-R": mt_ref,
        "-V": vcf,
        "-O": out_fn,
        "--autosomal-coverage": autosomal_coverage,
    }

    gatk_exec.run(subcommand="NuMTFilterTool", **params)


def _normalize_vcf(vcf: str, out_fn: str, mt_ref: str, gatk_exec: Executable) -> None:
    """Normalize variants calls using gatk LeftAlignAndTrimVariants."""

    params = {
        "-R": mt_ref,
        "-V": vcf,
        "-O": out_fn,
        "--split-multi-allelics": True,
        "--dont-trim-alleles": True,
        "--keep-original-ac": True,
    }

    gatk_exec.run(subcommand="LeftAlignAndTrimVariants", **params)


def _remove_non_pass(vcf: str, out_fn: str, gatk_exec: Executable) -> None:
    """Remove non pass variants using gatk SelectVariants."""

    params = {
        "-V": vcf,
        "-O": out_fn,
        "--exclude-filtered": True,
    }

    gatk_exec.run(subcommand="SelectVariants", **params)


def do_postprocess(
    vcf: str,
    stats: str = None,
    mt_ref: str = "rcrs",
    f_score_beta: float = 1.0,
    contamination_filter: bool = False,
    max_alt_allele_count: int = 4,
    vaf_treshold: float = 0,
    autosomal_coverage: float = 0,
    blacklisted_sites: str = None,
    remove_non_pass: bool = True,
    normalize: bool = True,
    out_dir: str = None,
    prefix: str = None,
    gatk_path: str = "gatk",
    haplocheck_path: str = "haplocheck",
    verbose: bool = False,
) -> dict:
    """Filter and normalize raw variant calls.

    Args:
        vcf (str): Path to VCF file
        stats (str, optional): Path to VCF stats file. Defaults to None.
        mt_ref (str, optional): Mitochondrial reference. Defaults to "rcrs".
        f_score_beta (float, optional): F score beta. Defaults to 1.
        contamination_filter (bool, optional): Activate contamination filter. Defaults to False.
        max_alt_allele_count (int, optional): Max alternate allele count. Defaults to 4.
        vaf_treshold (float, optional): VAF treshold. Defaults to 0.
        autosomal_coverage (float, optional): Median autosomal coverage. If set, NuMTs will be filtered. Defaults to 0.
        blacklisted_sites (str, optional): Path to BED file containing custom blacklisted sites. Defaults to None.
        remove_non_pass (bool, optional): Remove non passing variants. Defaults to True.
        normalize (bool, optional): Normalize variants. Defaults to True.
        out_dir (str, optional): Output directory. Defaults to None.
        prefix (str, optional): Prefix. Defaults to None.
        gatk_path (str, optional): Path to GATK executable. Defaults to "gatk".
        haplocheck_path (str, optional): Path to Haplocheck executable. Defaults to "haplocheck".
        verbose (bool, optional): Verbosity. Defaults to False.

    Returns:
        dict: Main output file paths
    """
    haplocheck = Executable(haplocheck_path, verbose)
    gatk = Executable(gatk_path, verbose)

    if not prefix:
        prefix = get_file_basename(vcf)

    if not out_dir:
        out_dir = get_file_directory(vcf)

    os.makedirs(out_dir, exist_ok=True)

    if not blacklisted_sites:
        blacklisted_sites = MT_BLACKLIST[mt_ref.lower()]

    mt_ref_fasta = MT_REFS[mt_ref.lower()]

    postprocessed_vcf = create_output_path(prefix, out_dir, "_postprocessed", ".vcf")

    # If contamination filter is enabled, estimate contamination, else set to 0.0
    contamination = (
        _get_contamination(vcf=vcf, out_dir=out_dir, haplocheck_exec=haplocheck)
        if contamination_filter
        else 0.0
    )

    # Check if stats file exists
    if not stats:
        stats = create_output_path(
            get_file_basename(vcf), get_file_directory(vcf), "", ".vcf.stats"
        )
        if not os.path.exists(stats):
            logging.error("VCF stats not found, please provide path to VCF stats.")

    # Initial filter
    logging.info("Filtering variants by parameters...")
    final_vcf = create_output_path(prefix, out_dir, "_filtered", ".vcf")
    _filter_by_params(
        vcf=vcf,
        mt_ref=mt_ref_fasta,
        stats=stats,
        vaf_treshold=vaf_treshold,
        max_alt_allele_count=max_alt_allele_count,
        f_score_beta=f_score_beta,
        contamination=contamination,
        out_fn=final_vcf,
        gatk_exec=gatk,
    )

    # Blacklist filtering
    logging.info("Filtering blacklisted sites...")
    blacklist_vcf = create_output_path(prefix, out_dir, "_blacklisted", ".vcf")

    _filter_blacklist(
        vcf=final_vcf,
        out_fn=blacklist_vcf,
        blacklisted_sites=blacklisted_sites,
        gatk_exec=gatk,
    )

    final_vcf = blacklist_vcf

    # Filtering of NuMTs
    if autosomal_coverage != 0:
        logging.info("Filtering Numts...")
        filtered_numt = create_output_path(prefix, out_dir, "_numt", ".vcf")
        _filter_numts(
            vcf=final_vcf,
            mt_ref=mt_ref_fasta,
            out_fn=filtered_numt,
            gatk_exec=gatk,
            autosomal_coverage=autosomal_coverage,
        )

        final_vcf = filtered_numt

    # Normalize
    if normalize:
        logging.info("Splitting multi-allelic sites and left-aligning variant calls...")
        normalized_vcf = create_output_path(prefix, out_dir, "_normalized", ".vcf")

        _normalize_vcf(
            vcf=final_vcf, out_fn=normalized_vcf, mt_ref=mt_ref_fasta, gatk_exec=gatk
        )

        final_vcf = normalized_vcf

    # Remove non pass variants
    if remove_non_pass:
        logging.info("Removing non pass variants...")
        pass_vcf = create_output_path(prefix, out_dir, "_pass", ".vcf")
        _remove_non_pass(vcf=final_vcf, out_fn=pass_vcf, gatk_exec=gatk)

        final_vcf = pass_vcf

    postprocessed_vcf = create_output_path(prefix, out_dir, "_postprocessed", ".vcf")
    postprocessed_vcf_idx = create_output_path(
        prefix, out_dir, "_postprocessed", ".vcf.idx"
    )

    os.rename(final_vcf, postprocessed_vcf)
    os.rename(f"{final_vcf}.idx", postprocessed_vcf_idx)

    # Collect outputs
    output_paths = {
        "postprocessed_vcf": postprocessed_vcf,
        "postprocessed_vcf_idx": postprocessed_vcf_idx,
    }

    # Check if output files exist
    if check_files_exist(list(output_paths.values())):
        logging.info(f"Variants postprocessing completed successfully.")
    else:
        logging.error("Some output files are missing! Please rerun the analysis.")
        sys.exit(1)

    return output_paths
