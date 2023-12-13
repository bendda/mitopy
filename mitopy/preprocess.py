from .utils import (
    check_bam_sorted,
    sort_bam,
    check_files_exist,
    index_bam,
    get_file_extension,
    create_sequence_dict,
    index_fasta,
    get_file_basename,
    get_file_directory,
    get_mt_contig_name,
    create_output_path,
)
import logging
import os
from .executable import Executable
import sys


def _subset_bam_chrm(
    bam: str,
    out_fn: str,
    gatk_exec: Executable,
    contig_name: str,
    reference_fa: str = None,
) -> None:
    """Subset BAM file to mt contig using gatk PrintReads."""

    params = {
        "-L": contig_name,
        "-I": bam,
        "-O": out_fn,
        "--read-filter": [
            "MateOnSameContigOrNoMappedMateReadFilter",
            "MateUnmappedAndUnmappedReadFilter",
        ],
    }

    if reference_fa:
        params["-R"] = reference_fa

    gatk_exec.run(subcommand="PrintReads", **params)


def _generate_ubam(bam: str, out_fn: str, gatk_exec: Executable) -> None:
    """Unalign mt reads using gatk RevertSam."""

    params = {
        "-I": bam,
        "-O": out_fn,
        "--VALIDATION_STRINGENCY": "LENIENT",
        "--ATTRIBUTE_TO_CLEAR": ["FT", "CO"],
        "--RESTORE_ORIGINAL_QUALITIES": False,
    }

    gatk_exec.run(subcommand="RevertSam", **params)


def do_preprocess(
    bam: str,
    bai: str = None,
    reference_fa: str = None,
    contig_name: str = None,
    out_dir: str = None,
    prefix: str = None,
    gatk_path: str = "gatk",
    verbose: bool = True,
) -> dict:
    """Preprocess input bam file for mitochondrial for mitochondrial variant calling and analysis.

    Args:
        bam (str): Path to input BAM file
        bai (str, optional): Path to BAI index file. Defaults to None.
        reference_fa (str, optional): Path to reference FASTA. Defaults to None.
        contig_name (str, optional): Name of the mitochondrial contig. Defaults to None.
        out_dir (str, optional): Output directory. Defaults to None.
        prefix (str, optional): Prefix. Defaults to None.
        gatk_path (str, optional): Path to GATK executable. Defaults to "gatk".
        verbose (bool, optional): Verbosity. If True, record logs of underlying tools. Defaults to True.

    Returns:
        dict: Main output file paths
    """
    gatk = Executable(gatk_path, verbose)

    ext = get_file_extension(bam)

    if not prefix:
        prefix = get_file_basename(bam)

    if not out_dir:
        out_dir = get_file_directory(bam)

    os.makedirs(out_dir, exist_ok=True)

    logging.info("Checking required input files...")

    # Check if input BAM is sorted
    if not check_bam_sorted(bam):
        logging.info("Input BAM file is not sorted. Sorting the BAM file...")
        sort_bam(bam)

    # Check if index exists
    index_exists = (
        bai
        or os.path.exists(f"{bam}.bai")
        or os.path.exists(f"{os.path.splitext(bam)[0]}.bai")
        or os.path.exists(f"{bam}.crai")
        or os.path.exists(f"{os.path.splitext(bam)[0]}.crai")
    )

    if not index_exists:
        logging.info(
            f"The BAM/CRAM index file not found/provided. Creating BAI/CRAI index file using pysam..."
        )
        index_bam(bam)

    # If CRAM, check if reference FASTA is provided (and dict, fai files exist)
    if ext == ".cram":
        if not reference_fa:
            logging.error(
                "Reference FASTA is required when input is a CRAM file. Please provide reference FASTA."
            )
            sys.exit(1)

        if not (os.path.exists(f"{os.path.splitext(reference_fa)[0]}.dict")):
            logging.info(
                "The reference dictionary file not found. Creating reference dictionary file..."
            )
            create_sequence_dict(reference_fa)

        if not (os.path.exists(f"{reference_fa}.fai")):
            logging.info(
                f"The reference FASTA index file not found. Creating FAI index file..."
            )
            index_fasta(reference_fa)

    # Detect contig name if not provided
    if not contig_name:
        try:
            logging.info(
                "Contig name not provided. Detecting mitochondrial contig name from input BAM/CRAM file..."
            )
            contig_name = get_mt_contig_name(bam)
        except ValueError:
            logging.error(
                "Mitochondrial contig name could not be retrieved. It is expected to be either MT or chrM based on used human reference."
            )
            sys.exit(1)

    # Subset input BAM to mitochondrial contig
    logging.info(
        f"Subsetting bam file to keep only reads mapped to {contig_name} contig..."
    )
    subset_out = create_output_path(prefix, out_dir, "_chrM", ".bam")

    _subset_bam_chrm(
        bam=bam,
        out_fn=subset_out,
        gatk_exec=gatk,
        contig_name=contig_name,
        reference_fa=reference_fa,
    )

    # Generate unmapped BAM
    logging.info("Generating unmapped BAM file...")
    revert_out = create_output_path(prefix, out_dir, "_unmapped", ".bam")

    _generate_ubam(bam=subset_out, out_fn=revert_out, gatk_exec=gatk)

    # Collect outputs
    output_paths = {"unmapped_bam": revert_out}

    # Check if output files exist
    if check_files_exist(list(output_paths.values())):
        logging.info(f"Preprocessing of {bam} completed successfully.")
    else:
        logging.error("Some output files are missing! Please rerun the analysis.")
        sys.exit(1)

    return output_paths
