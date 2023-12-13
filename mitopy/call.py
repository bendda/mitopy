from .executable import Executable
from .utils import (
    get_file_basename,
    get_file_directory,
    create_output_path,
    check_files_exist,
)
import logging
import os
import sys
from .constants import MT_REFS

NON_CONTROL_REGION = "chrM:576-16024"
CONTROL_REGION = "chrM:8025-9144"


def do_call(
    bam: str,
    mt_ref: str = "rcrs",
    m2_extra_args: str = "",
    out_dir: str = None,
    prefix: str = None,
    gatk_path: str = "gatk",
    shifted: bool = False,
    verbose: bool = False,
) -> dict:
    """Call mitochondrial variants using Mutect2.

    Args:
        bam (str): Path to BAM file
        mt_ref (str, optional): Mitochondrial reference. Defaults to "rcrs".
        m2_extra_args (str, optional): Extra args to pass onto Mutect2. Defaults to "".
        out_dir (str, optional): Output directory. Defaults to None.
        prefix (str, optional): Prefix. Defaults to None.
        gatk_path (str, optional): Path to GATK executable. Defaults to "gatk".
        shifted (bool, optional): Shifted mode. If True, variants are called against shifted mitochondrial reference. Defaults to False.
        verbose (bool, optional): Verbosity. If True, record logs of underlying tools. Defaults to False.

    Returns:
        dict: Main output file paths
    """
    gatk = Executable(gatk_path, verbose)

    if not prefix:
        prefix = get_file_basename(bam)

    if not out_dir:
        out_dir = get_file_directory(bam)

    os.makedirs(out_dir, exist_ok=True)

    if shifted:
        prefix += "_shifted"

    mt_reference_fasta = (
        MT_REFS[mt_ref.lower()] if not shifted else MT_REFS[f"{mt_ref.lower()}_shifted"]
    )

    # Call variants
    output_fn = create_output_path(prefix, out_dir, "", ".vcf")

    params = {
        "-I": bam,
        "-R": mt_reference_fasta,
        "-O": output_fn,
        "--read-filter": [
            "MateOnSameContigOrNoMappedMateReadFilter",
            "MateUnmappedAndUnmappedReadFilter",
        ],
        "--annotation": "StrandBiasBySample",
        "--mitochondria-mode": True,
        "--max-reads-per-alignment-start": 75,
        "--max-mnp-distance": 0,
        "-L": CONTROL_REGION if shifted else NON_CONTROL_REGION,
    }

    logging.info("Calling variants with Mutect2 in mitochondria mode...")
    gatk.run(m2_extra_args, subcommand="Mutect2", **params)

    # Collect outputs
    output_paths = {
        "raw_vcf": output_fn,
        "raw_vcf_index": f"{output_fn}.idx",
        "raw_vcf_stats": f"{output_fn}.stats",
    }

    # Check if output files exist
    if check_files_exist(list(output_paths.values())):
        logging.info(f"Variant calling completed successfully.")
    else:
        logging.error("Some output files are missing! Please rerun the analysis.")
        sys.exit(1)

    return output_paths
