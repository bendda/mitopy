import os
import logging
import pysam
from pathlib import Path
import subprocess


def check_files_exist(files: list | str, verbose: bool = False) -> bool:
    """Check if the files exist."""

    files_list = files if isinstance(files, list) else [files]
    for f in files_list:
        if not os.path.isfile(f):
            if verbose:
                logging.info(f"File {f} not found.")
            return False
    return True


def check_bam_sorted(input_bam: str) -> bool:
    """Check if BAM file is sorted."""

    # Workaround to not print missing index file warning
    pysam.set_verbosity(0)

    with pysam.AlignmentFile(input_bam, "rb") as bam:
        return bam.header["HD"]["SO"] == "coordinate"


def sort_bam(input_bam: str) -> str:
    """Sort BAM file."""

    basename = os.path.splitext(input_bam)[0]
    ext = os.path.splitext(input_bam)[1]
    output_fn = f"{basename}.sorted{ext}"

    pysam.sort(input_bam, "-o", output_fn)
    return output_fn


def index_bam(input_bam: str) -> str:
    """Index BAM/CRAM file."""

    ext = get_file_extension(input_bam)

    pysam.index(input_bam)
    return f"{input_bam}.bai" if ext == ".bam" else f"{input_bam}.crai"


def create_sequence_dict(input_fasta: str) -> str:
    """Create sequence dictionary for FASTA file."""

    basename = os.path.splitext(input_fasta)[0]
    cmd = f"gatk CreateSequenceDictionary -R {input_fasta}"
    subprocess.run(cmd, shell=True)
    return f"{basename}.dict"


def index_fasta(input_fasta: str) -> str:
    """Index FASTA file."""

    pysam.faidx(input_fasta)
    return f"{input_fasta}.fai"


def get_mt_contig_name(input_bam: str) -> str:
    """Get mitochondrial contig name from BAM file."""

    with pysam.AlignmentFile(input_bam, "rb") as bam:
        chroms = [str(record.get("SN")) for record in bam.header["SQ"]]
        mito_contig = {"MT", "chrM"}.intersection(chroms)
        if mito_contig is not None:
            return "".join(mito_contig)
        else:
            raise ValueError(
                "Mitochondrial contig (MT or chrM) not found in the BAM file."
            )


def get_file_directory(file_path: str) -> str:
    """Get directory of the file."""
    return str(Path(file_path).parent)


def get_file_basename(file_path: str) -> str:
    """Get basename of the file."""
    return str(Path(file_path).stem)


def get_file_extension(file_path: str) -> str:
    """Get file extension."""
    return str(Path(file_path).suffix)


def create_output_path(prefix, out_dir, suffix, ext):
    """Create file path."""
    return f"{out_dir}/{prefix}{suffix}{ext}"
