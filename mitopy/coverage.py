import pandas as pd
import plotly.graph_objects as go
import os
from .utils import (
    get_file_basename,
    get_file_directory,
    create_output_path,
    check_files_exist,
)
from .executable import Executable
import sys
import logging

pd.options.mode.chained_assignment = None


def _get_mosdepth_pb_coverage(
    bam: str, prefix: str, mosdepth_exec: Executable
) -> pd.DataFrame:
    """Compute per-base coverage using mosdepth."""

    mosdepth_exec.run(prefix, bam)

    # Load to dataframe
    pb_cov = pd.read_csv(
        f"{prefix}.per-base.bed.gz",
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "coverage"],
    )

    # Convert to VCF coordinate system
    # 0-based, half-open [start-1, end) --> 1-based, closed [start, end]
    pb_cov["start"] = pb_cov["start"] + 1

    return pb_cov


def _shift_back(pos: int) -> int:
    """Shift the position back to canonical coordinates."""
    return pos + 8000 if pos < 8570 else pos - 8569


def _plot_coverage(coverage_csv: pd.DataFrame) -> go.Figure:
    """Create coverage plot."""

    fig = go.Figure()

    fig.add_trace(
        go.Scatter(
            x=coverage_csv["start"],
            y=coverage_csv["coverage"],
            mode="lines",
            name="Coverage per Base",
        )
    )

    # Set axis labels and title
    fig.update_layout(
        xaxis=dict(title="Genomic Position"),
        yaxis=dict(title="Coverage"),
        title="Coverage per Base",
    )

    return fig


def do_coverage(
    mt_bam: str,
    shifted_mt_bam: str,
    mt_bai: str = None,
    shifted_mt_bai: str = None,
    out_dir: str = None,
    prefix: str = None,
    create_plot: bool = True,
    mosdepth_path: str = "mosdepth",
    verbose: bool = False,
) -> dict:
    """Calculate and combine per-base coverage of control and non-control region.

    Args:
        mt_bam (str): Path to BAM file (canonical alignment)
        shifted_mt_bam (str): Path to BAM file (shifted alignment)
        mt_bai (str, optional): BAM index file. Defaults to None.
        shifted_mt_bai (str, optional): Shifted BAM index file. Defaults to None.
        out_dir (str, optional): Output directory. Defaults to None.
        prefix (str, optional): Prefix. Defaults to None.
        create_plot (bool, optional): Create coverage plot. Defaults to True.
        mosdepth_path (str, optional): Path to mosdepth executable. Defaults to "mosdepth".
        verbose (bool, optional): Verbosity. Defaults to False.

    Returns:
        dict: Main output file paths
    """
    mosdepth = Executable(mosdepth_path, verbose)

    if not prefix:
        prefix = get_file_basename(mt_bam)

    if not out_dir:
        out_dir = get_file_directory(mt_bam)

    os.makedirs(out_dir, exist_ok=True)

    mosdepth_dir = os.path.join(out_dir, "tmp")
    os.makedirs(mosdepth_dir, exist_ok=True)

    # Get per-base coverage for non-control region
    logging.info(f"Getting per base coverage for non-control region using mosdepth...")
    df = _get_mosdepth_pb_coverage(
        mt_bam, f"{mosdepth_dir}/non_control", mosdepth_exec=mosdepth
    )

    # Get per-base coverage for control region
    logging.info(f"Getting per base coverage for control region using mosdepth...")
    df_shifted = _get_mosdepth_pb_coverage(
        shifted_mt_bam, f"{mosdepth_dir}/control", mosdepth_exec=mosdepth
    )

    logging.info(f"Combining per base coverage from control and non-control region...")

    # Subset non-control region
    non_control = df[(df["start"] > 576) & (df["end"] < 16024)]

    # Subset control region and shift-back to canonical coordinates
    control = df_shifted[(df_shifted["start"] > 8023) & (df_shifted["end"] < 9147)]
    control["start"] = control["start"].apply(_shift_back)
    control["end"] = control["end"].apply(_shift_back)

    start = control[control["start"] < 8000]
    end = control[control["start"] > 8000]

    # Combine per-base coverage from control and non-control regions
    coverage_csv = create_output_path(prefix, out_dir, "_coverage", ".csv")
    combined_perbase = pd.concat([start, non_control, end], ignore_index=True)
    combined_perbase.to_csv(coverage_csv, index=False)

    output_paths = {
        "coverage_csv": coverage_csv,
    }

    if create_plot:
        logging.info("Creating per-base coverage plot...")
        fig = _plot_coverage(combined_perbase)

        # Save as interactive (html)
        coverage_html = create_output_path(prefix, out_dir, "_coverage", ".html")
        fig.write_html(coverage_html)
        output_paths["coverage_html"] = coverage_html

    # Check if output files exist
    if check_files_exist(list(output_paths.values())):
        logging.info(f"Calculating combined coverage per-base completed successfully.")
    else:
        logging.error("Some output files are missing! Please rerun the analysis.")
        sys.exit(1)

    return output_paths
