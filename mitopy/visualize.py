import pandas as pd
import plotly.graph_objects as go
from pysam import VariantFile
import sys
import os
from .utils import (
    get_file_directory,
    create_output_path,
    get_file_basename,
    check_files_exist,
)
from .constants import VIS_RESOURCES
import logging


color_type_mapping = {
    "tRNA": "green",
    "rRNA": "yellow",
    "REG": "black",
    "CDS": "darkred",
}

MT_LENGTH = 16569


def _convert_to_polar(pos: int) -> float:
    """Convert position to polar coordinate system."""
    return (pos / MT_LENGTH) * 360


def _parse_vcf(vcf: str) -> dict:
    """Parse VCF file"""
    vcf_file = VariantFile(vcf)
    sample_name = vcf_file.header.samples[0]

    var_pos = []
    var_af = []
    var_ref = []
    var_alt = []

    for rec in vcf_file.fetch():
        var_pos.append(rec.pos)
        var_af.append(rec.samples[0]["AF"][0])
        var_ref.append(rec.ref)
        var_alt.append(rec.alts)

    return {
        "position": var_pos,
        "vaf": var_af,
        "sample": sample_name,
        "ref": var_ref,
        "alt": var_alt,
    }


def _create_feature_trace(mito_map: list, feature_type: str) -> go.Figure:
    """Create mito feature trace."""

    feature = [tup for tup in mito_map if tup[5] == feature_type]

    return go.Barpolar(
        r=[tup[0] for tup in feature],
        theta=[tup[1] for tup in feature],
        width=[tup[2] for tup in feature],
        base=[tup[3] for tup in feature],
        marker_color=[tup[4] for tup in feature],
        text=[tup[6] for tup in feature],
        hoverinfo="text",
        name=feature_type,
        showlegend=True,
    )


def _plot_mtbase(split: bool = False) -> go.Figure:
    """Plot mitochondrial genome map."""

    mito_annot = pd.read_csv(VIS_RESOURCES["mitomap"])
    num_features = mito_annot.shape[0]
    trace_width = 30
    trace_start = 70

    fig = go.Figure()

    widths = list(map(_convert_to_polar, mito_annot["length"]))
    feature_types = mito_annot["type"].tolist()
    colors = [color_type_mapping[type] for type in feature_types]
    names = mito_annot["gene_name"].tolist()

    radius = [
        trace_width if not split else trace_width / 2 for _ in range(num_features)
    ]
    base = [
        trace_start if (strand == "+" or not split) else (trace_start + trace_width / 2)
        for strand in mito_annot["strand"]
    ]

    # Calculate thetas (positions of angular axis)
    theta = [0.5 * w for w in widths]
    for index, w in enumerate(widths):
        for subsequent_index in range(index + 1, len(widths)):
            theta[subsequent_index] += w

    # MT coordinates
    mt_coordinates = list(range(0, MT_LENGTH + 1, 1000))
    mt_coordinates_polar = list(map(_convert_to_polar, mt_coordinates))
    mt_coordinates_labels = [f"{c} bp" for c in mt_coordinates]

    mito_map = list(zip(radius, theta, widths, base, colors, feature_types, names))

    # Add border trace
    border_trace = go.Barpolar(
        r=[trace_width, 0.3, 0.3],
        theta=[0, 0, 0],
        width=[360, 360, 360],
        base=[trace_start, trace_start, trace_start + trace_width],
        marker_color=["#EFEFEF", "black", "black"],
        hoverinfo="skip",
        showlegend=False,
    )

    fig.add_trace(border_trace)

    # Add feature traces
    for feature_type in color_type_mapping.keys():
        feature_trace = _create_feature_trace(mito_map, feature_type)
        fig.add_trace(feature_trace)

    # Create layout
    fig.update_layout(
        title="Mitochondria rCRS",
        width=1000,
        height=1000,
        hovermode="closest",
        polar=dict(
            radialaxis=dict(
                showgrid=False, showline=False, showticklabels=False, range=[0, 100]
            ),
            angularaxis=dict(
                dtick=20,
                tick0=0,
                showgrid=False,
                showline=False,
                linecolor="black",
                direction="clockwise",
                tickvals=mt_coordinates_polar,
                ticktext=mt_coordinates_labels,
                tickangle=0,
                ticks="outside",
            ),
            bgcolor="white",
            hole=0.3,
        ),
    )

    return fig


def _add_variant_trace(fig: go.Figure, vcf: str) -> go.Figure:
    """Add variants to the figure."""

    vcf_info = _parse_vcf(vcf)

    trace_width = 20
    trace_start = 40

    # Scale
    scaled_var_af = [trace_start + vaf * trace_width for vaf in vcf_info["vaf"]]
    var_pos_polar = list(map(_convert_to_polar, vcf_info["position"]))

    annotation = list(
        zip(vcf_info["position"], vcf_info["ref"], vcf_info["alt"], vcf_info["vaf"])
    )

    # Define variant trace
    var_trace = go.Scatterpolar(
        r=scaled_var_af,
        theta=var_pos_polar,
        mode="markers",
        marker=dict(color="black"),
        customdata=annotation,
        hovertemplate="Position: %{customdata[0]}<br>"
        + "Reference Allele: %{customdata[1]}<br>"
        + "Alternate Allele: %{customdata[2]}<br>"
        + "Heteroplasmy Fraction: %{customdata[3]:.2f}<extra></extra>",
        showlegend=True,
        name="Variants",
    )

    # Define border trace
    border_trace = go.Barpolar(
        r=[trace_width, 0.3, 0.3],
        theta=[0, 0, 0],
        width=[360, 360, 360],
        base=[trace_start, trace_start, trace_start + trace_width],
        marker_color=["#EFEFEF", "black", "black"],
        hoverinfo="skip",
        showlegend=False,
    )

    # Add traces
    fig.add_trace(border_trace)
    fig.add_trace(var_trace)

    # Add sample name
    fig.update_layout(title=vcf_info["sample"])

    return fig


def _add_coverage_trace(fig: go.Figure, coverage_csv: str) -> go.Figure:
    """Add coverage plot to figure."""

    trace_width = 20
    trace_start = 0

    # Load coverage CSV
    coverage_info = pd.read_csv(coverage_csv)
    coverage = coverage_info["coverage"]
    position_bp = coverage_info["start"]

    annotation = list(zip(position_bp, coverage))

    # Scale to polar coordinates
    max_cov = max(coverage)
    scaled_coverage = [(cov * trace_width) / max_cov for cov in coverage]
    position_polar = list(map(_convert_to_polar, position_bp))

    # Define coverage trace
    coverage_trace = go.Scatterpolar(
        theta=position_polar,
        r=scaled_coverage,
        mode="lines",
        name="Coverage per Base",
        customdata=annotation,
        line_color="#E3735E",
        hoveron="points",
        fill="toself",
        hovertemplate="Position: %{customdata[0]}<br>Coverage: %{customdata[1]:.2f}<extra></extra>",
    )

    # Define border trace
    border_trace = go.Barpolar(
        r=[0.3, 0.3],
        theta=[0, 0],
        width=[360, 360],
        base=[trace_start, trace_start + trace_width],
        marker_color=["black", "black"],
        hoverinfo="skip",
        showlegend=False,
    )

    # Add traces
    fig.add_trace(border_trace)
    fig.add_trace(coverage_trace)

    return fig


def do_visualize(
    vcf: str,
    coverage_csv: str = None,
    split_strands: bool = False,
    save_as_png: bool = False,
    out_dir: str = None,
    prefix: str = None,
) -> dict:
    """Visualize variant calls.

    Args:
        vcf (str): Path to input VCF file
        coverage_csv (str, optional): Path to CSV containing coverage. Defaults to None.
        split_strands (bool, optional): Split strands on mito genome. Defaults to False.
        save_as_png (bool, optional): Save plot as PNG. Defaults to False.
        out_dir (str, optional): Output directory. Defaults to None.
        prefix (str, optional): Prefix. Defaults to None.

    Returns:
        dict: Main output file paths
    """

    if not out_dir:
        out_dir = get_file_directory(vcf)

    if not prefix:
        prefix = get_file_basename(vcf)

    os.makedirs(out_dir, exist_ok=True)

    # Plot mitochondrial genome map
    logging.info("Plotting mtiochondrial genome base...")
    fig = _plot_mtbase(split_strands)

    # Plot variants
    logging.info("Plotting variants...")
    fig = _add_variant_trace(fig, vcf)

    # Add coverage plot
    if coverage_csv:
        logging.info("Plotting per base coverage...")
        fig = _add_coverage_trace(fig, coverage_csv)

    # Save as HTML
    out_html = create_output_path(
        prefix=prefix, out_dir=out_dir, suffix="", ext=".html"
    )
    fig.write_html(out_html)

    output_paths = {
        "vis_html": out_html,
    }

    # Save as PNG
    if save_as_png:
        logging.info("Saving as PNG image...")
        out_png = create_output_path(
            prefix=prefix, out_dir=out_dir, suffix="", ext=".png"
        )
        fig.write_image(out_png)
        output_paths["vis_png"] = out_png

    # Check if output files exist
    if check_files_exist(list(output_paths.values())):
        logging.info(f"Variant visualization completed successfully.")
    else:
        logging.error("Some output files are missing! Please rerun the analysis.")
        sys.exit(1)

    return output_paths
