from .executable import Executable
from .utils import (
    get_file_basename,
    get_file_directory,
    create_output_path,
    check_files_exist,
)
import os
import logging
import sys


PHYLOTREES = {"rcrs": "phylotree-rcrs@17.2", "rsrs": "phylotree-rsrs@17.1"}


def do_identify_haplogroup(
    vcf: str,
    mt_ref: str = "rcrs",
    haplogrep3_path: str = "haplogrep3",
    prefix: str = None,
    out_dir: str = None,
    verbose: bool = False,
) -> dict:
    """Identify haplogroup using Haplogrep3

    Args:
        vcf (str): path to VCF file
        mt_ref (str, optional): Mitochondrial reference. Defaults to "rcrs".
        haplogrep3_path (str, optional): Haplogrep3 path. Defaults to "haplogrep3".
        prefix (str, optional): Prefix. Defaults to None.
        out_dir (str, optional): Output directory. Defaults to None.
        verbose (bool, optional): Verbosity. Defaults to False.

    Returns:
        dict: Main output file names
    """

    haplogrep3 = Executable(haplogrep3_path, verbose)

    if not prefix:
        prefix = get_file_basename(vcf)

    if not out_dir:
        out_dir = get_file_directory(vcf)

    os.makedirs(out_dir, exist_ok=True)

    # Select phylotree
    tree = PHYLOTREES[mt_ref.lower()]

    # Download phylotree for rsrs
    if mt_ref.lower() == "rsrs":
        logging.info("Downloading phylotree for RSRS reference...")
        haplogrep3.run(tree, subcommand="install-tree")

    # Classify
    logging.info("Identifying haplogroup using Haplogrep3...")
    haplo_out = create_output_path(prefix, out_dir, "_haplogroup", ".txt")
    params = {"--tree": tree, "--in": vcf, "--out": haplo_out}
    haplogrep3.run(subcommand="classify", **params)

    output_paths = {
        "haplogroups": haplo_out,
    }
    # Check if output files exist
    if check_files_exist(list(output_paths.values())):
        logging.info(f"Haplogroup identification completed successfully.")
    else:
        logging.error("Some output files are missing! Please rerun the analysis.")
        sys.exit(1)

    return output_paths
