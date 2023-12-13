import logging
from .executable import Executable
from .utils import (
    get_file_basename,
    get_file_directory,
    create_output_path,
    check_files_exist,
)
import os
from .constants import MT_REFS
import sys


def _liftover_shifted(
    vcf_shifted: str,
    mt_ref: str,
    shift_back_chain: str,
    out_fn: str,
    rejected_out: str,
    gatk_exec: Executable,
) -> None:
    """Liftover VCF file using gatk LiftoverVcf."""

    params = {
        "-I": vcf_shifted,
        "-R": mt_ref,
        "-O": out_fn,
        "-C": shift_back_chain,
        "--REJECT": rejected_out,
    }

    gatk_exec.run(subcommand="LiftoverVcf", **params)


def _merge_vcfs(
    vcf: str, vcf_shifted_back: str, out_fn: str, gatk_exec: Executable
) -> None:
    """Merge VCF files using gatk MergeVcfs."""

    params = {"-I": [vcf, vcf_shifted_back], "-O": out_fn}

    gatk_exec.run(subcommand="MergeVcfs", **params)


def _merge_stats(
    stats: str, stats_shifted: str, out_fn: str, gatk_exec: Executable
) -> None:
    """Merge VCF stats using gatk MergeMutectStats."""

    params = {"-stats": [stats, stats_shifted], "-O": out_fn}

    gatk_exec.run(subcommand="MergeMutectStats", **params)


def do_merge(
    vcf: str,
    vcf_shifted: str,
    mt_ref: str = "rcrs",
    stats: str = None,
    stats_shifted: str = None,
    out_dir: str = None,
    prefix: str = None,
    gatk_path: str = "gatk",
    verbose: bool = False,
) -> dict:
    """Merge variant calls and stats from control and non-control mt regions.

    Args:
        vcf (str): Path to VCF file (non-control region)
        vcf_shifted (str): Path to shifted VCF file (control region)
        mt_ref (str, optional): Mitochondrial reference. Defaults to "rcrs".
        stats (str, optional): VCF stats file. Defaults to None.
        stats_shifted (str, optional): Shifted VCF stats file. Defaults to None.
        out_dir (str, optional): Output directory. Defaults to None.
        prefix (str, optional): Prefix. Defaults to None.
        gatk_path (str, optional): Path to GATK executable. Defaults to "gatk".
        verbose (bool, optional): Verbosity. If true, record logs of underlying tools. Defaults to False.

    Returns:
        dict: Main output file paths
    """
    gatk = Executable(gatk_path, verbose)

    if not prefix:
        prefix = get_file_basename(vcf)

    if not out_dir:
        out_dir = get_file_directory(vcf)

    os.makedirs(out_dir, exist_ok=True)

    mt_ref_fasta = MT_REFS[mt_ref.lower()]

    # If stats files are not explicitly provided, check if they exist in directory of VCF
    if not stats:
        stats = create_output_path(
            get_file_basename(vcf), get_file_directory(vcf), "", ".vcf.stats"
        )
        if not os.path.exists(stats):
            logging.error("VCF stats not found, please provide path to vcf stats.")
            sys.exit(1)

    if not stats_shifted:
        stats_shifted = create_output_path(
            get_file_basename(vcf_shifted),
            get_file_directory(vcf_shifted),
            "",
            ".vcf.stats",
        )
        if not os.path.exists(stats_shifted):
            logging.error(
                "VCF stats not found, please provide path to shifted vcf stats."
            )
            sys.exit(1)

    shift_back_chain = MT_REFS[f"{mt_ref.lower()}_shift_back_chain"]

    # Shift back shifted VCF
    logging.info("Lifting over shifted VCF...")
    vcf_shifted_back = create_output_path(prefix, out_dir, "_shifted_back", ".vcf")
    vcf_rejected = create_output_path(prefix, out_dir, "_rejected", ".vcf")

    _liftover_shifted(
        vcf_shifted=vcf_shifted,
        mt_ref=mt_ref_fasta,
        shift_back_chain=shift_back_chain,
        out_fn=vcf_shifted_back,
        rejected_out=vcf_rejected,
        gatk_exec=gatk,
    )

    # Merge VCFs
    logging.info("Merging VCFs...")
    merged_vcf = create_output_path(prefix, out_dir, "_merged", ".vcf")

    _merge_vcfs(
        vcf=vcf, vcf_shifted_back=vcf_shifted_back, out_fn=merged_vcf, gatk_exec=gatk
    )

    # Merge VCF stats
    logging.info("Merging Mutect2 stats...")
    merged_stats = create_output_path(prefix, out_dir, "_merged", ".vcf.stats")
    _merge_stats(
        stats=stats, stats_shifted=stats_shifted, out_fn=merged_stats, gatk_exec=gatk
    )

    # Collect outputs
    output_paths = {
        "merged_vcf": merged_vcf,
        "merged_vcf_index": f"{merged_vcf}.idx",
        "merged_vcf_stats": f"{merged_vcf}.stats",
    }

    # Check if output files exist
    if check_files_exist(list(output_paths.values())):
        logging.info(f"Merging VCF files completed successfully.")
    else:
        logging.error("Some output files are missing! Please rerun the analysis.")
        sys.exit(1)

    return output_paths
