import logging

from .executable import Executable
from .utils import (
    get_file_basename,
    get_file_directory,
    create_output_path,
    check_files_exist,
)
from .constants import MT_REFS
import os
import sys


def _ubam_to_fastq(ubam: str, out_fn: str, gatk_exec: Executable) -> None:
    """Convert uBAM to FASTQ format using gatk SamToFastq."""

    params = {
        "-I": ubam,
        "-F": out_fn,
        "--INTERLEAVE": True,
        "--INCLUDE_NON_PF_READS": True,
    }

    gatk_exec.run(subcommand="SamToFastq", **params)


def _merge_alignment(
    bam: str, ubam: str, mt_ref: str, out_fn: str, gatk_exec: Executable
) -> None:
    """Merge alignment with uBAM using gatk MergeBamAlignment."""

    params = {
        "--ALIGNED_BAM": bam,
        "--UNMAPPED_BAM": ubam,
        "--REFERENCE_SEQUENCE": mt_ref,
        "--OUTPUT": out_fn,
        "--VALIDATION_STRINGENCY": "SILENT",
        "--EXPECTED_ORIENTATIONS": "FR",
        "--ATTRIBUTES_TO_RETAIN": "X0",
        "--ATTRIBUTES_TO_REMOVE": ["NM", "MD"],
        "--SORT_ORDER": "queryname",
        "--CLIP_ADAPTERS": False,
        "--MAX_RECORDS_IN_RAM": 2000000,
        "--MAX_INSERTIONS_OR_DELETIONS": -1,
        "--PRIMARY_ALIGNMENT_STRATEGY": "MostDistant",
        "--UNMAPPED_READ_STRATEGY": "COPY_TO_TAG",
        "--ALIGNER_PROPER_PAIR_FLAGS": True,
        "--UNMAP_CONTAMINANT_READS": True,
        "--ADD_PG_TAG_TO_READS": False,
    }

    gatk_exec.run(subcommand="MergeBamAlignment", **params)


def _mark_dup_sort(
    bam: str, out_fn: str, out_metrics: str, ncores: int, gatk_exec: Executable
) -> None:
    """Mark duplicates and coordinate-sort using gatk MarkDuplicatesSpark."""

    ncores = "*" if ncores == -1 else ncores
    params = {
        "-I": bam,
        "-O": out_fn,
        "-M": out_metrics,
        "--optical-duplicate-pixel-distance": 2500,
        "--create-output-bam-splitting-index": False,
        "--spark-runner": "LOCAL",
        "--spark-master": f"local[{ncores}]",
    }

    gatk_exec.run(subcommand="MarkDuplicatesSpark", **params)


def _mt_bwa_align(
    mt_fastq: str, mt_ref: str, out_fn: str, ncores: int, bwa_exec: Executable
) -> None:
    """Align reads to mitochondrial reference using bwa-mem2."""
    params = {
        "-p": True,
        "-v": 3,
        "-t": ncores,
        "-K": 100000000,
        "-Y": True,
        "-o": out_fn,
    }

    bwa_exec.run(mt_ref, mt_fastq, subcommand="mem", **params)


def do_align(
    ubam: str,
    mt_ref: str = "rcrs",
    out_dir: str = None,
    prefix: str = None,
    shifted: bool = False,
    ncores: int = 1,
    verbose: bool = False,
    gatk_path: str = "gatk",
    bwamem2_path: str = "bwa-mem2",
) -> dict:
    """Align unmapped BAM file containing mitochondrial reads to mitochondrial reference.

    Args:
        ubam (str): Path to uBAM
        mt_ref (str, optional): Mitochondrial reference. Defaults to "rcrs".
        out_dir (str, optional): Output directory. Defaults to None.
        prefix (str, optional): Prefix. Defaults to None.
        shifted (bool, optional): Shifted mode. If True, align against shifted mitochondrial reference. Defaults to False.
        ncores (int, optional): Number of cores. Defaults to 1.
        verbose (bool, optional): Verbosity. If True, record logs of underlying tools. Defaults to False.
        gatk_path (str, optional): Path to GATK executable. Defaults to "gatk".
        bwamem2_path (str, optional): Path to bwa-mem2 executable. Defaults to "bwa-mem2".

    Returns:
        dict: Main output file paths
    """

    gatk = Executable(gatk_path, verbose)
    bwamem2 = Executable(bwamem2_path, verbose)

    if not prefix:
        prefix = get_file_basename(ubam)

    if not out_dir:
        out_dir = get_file_directory(ubam)

    os.makedirs(out_dir, exist_ok=True)

    if shifted:
        prefix += "_shifted"

    mt_ref_fasta = (
        MT_REFS[mt_ref.lower()] if not shifted else MT_REFS[f"{mt_ref.lower()}_shifted"]
    )

    # Convert uBAM to FASTQ
    logging.info("Converting uBAM to FASTQ format for alignment...")
    mt_fastq = create_output_path(prefix=prefix, out_dir=out_dir, suffix="", ext=".fq")
    _ubam_to_fastq(ubam=ubam, out_fn=mt_fastq, gatk_exec=gatk)

    # Align to mitochondrial reference
    logging.info(f"Aligning to mitochondrial reference genome {mt_ref}...")
    aligned_sam = create_output_path(
        prefix=prefix, out_dir=out_dir, suffix="", ext=".sam"
    )

    _mt_bwa_align(
        mt_fastq=mt_fastq,
        mt_ref=mt_ref_fasta,
        out_fn=aligned_sam,
        ncores=ncores,
        bwa_exec=bwamem2,
    )

    # Merge alignment with uBAM
    logging.info("Merging aligned SAM with uBAM...")
    merged_out = create_output_path(prefix, out_dir, "_merged", ".bam")
    _merge_alignment(
        bam=aligned_sam,
        ubam=ubam,
        mt_ref=mt_ref_fasta,
        out_fn=merged_out,
        gatk_exec=gatk,
    )

    # Mark duplicates and coordinate-sort
    logging.info("Marking duplicates and coordinate-sorting...")
    md_out = create_output_path(prefix, out_dir, "_dedup", ".bam")
    md_metrics = create_output_path(prefix, out_dir, "", ".dedup.metrics.txt")
    _mark_dup_sort(
        bam=merged_out,
        out_fn=md_out,
        out_metrics=md_metrics,
        ncores=ncores,
        gatk_exec=gatk,
    )

    # Collect main outputs
    output_paths = {
        "dedup_sorted_bam": md_out,
        "dedup_sorted_bai": f"{md_out}.bai",
    }

    # Check if output files exist
    if check_files_exist(list(output_paths.values())):
        logging.info(
            f"Alignment of {ubam} to {'shifted' if shifted else ''} mitochondrial reference {mt_ref} completed successfully."
        )
    else:
        logging.error("Some output files are missing! Please rerun the analysis.")
        sys.exit(1)

    return output_paths
