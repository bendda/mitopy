from .preprocess import do_preprocess
from .align import do_align
from .call import do_call
from .merge import do_merge
from .postprocess import do_postprocess
from .annotate import do_annotate
from .visualize import do_visualize
from .haplogroup import do_identify_haplogroup
from .coverage import do_coverage
import logging
import os
from .utils import get_file_basename, get_file_directory, check_files_exist
import shutil
import sys


def do_run_pipeline(
    bam: str,
    bai: str = None,
    mt_ref: str = "rcrs",
    reference_fa: str = None,
    contig_name: str = None,
    out_dir: str = None,
    tmp_dir: str = None,
    remove_tmp: bool = False,
    prefix: str = None,
    ncores: int = 1,
    m2_extra_args: str = None,
    f_score_beta: float = 1,
    contamination_filter: bool = True,
    max_alt_allele_count: int = 4,
    vaf_treshold: float = 0,
    autosomal_coverage: float = 0,
    blacklisted_sites: str = None,
    remove_non_pass: bool = True,
    normalize: bool = True,
    min_hom_treshold: float = 0.95,
    split_strands: bool = True,
    population_freqs: bool = True,
    patho_predictions: bool = True,
    phenotype_annot: bool = True,
    conservation_scores: bool = True,
    save_as_png: bool = True,
    create_annotation_report: bool = True,
    snpeff_path: str = "snpeff",
    snpsift_path: str = "snpsift",
    haplogrep3_path: str = "haplogrep3",
    gatk_path: str = "gatk",
    bwamem2_path: str = "bwa-mem2",
    haplocheck_path: str = "haplocheck",
    mosdepth_path: str = "mosdepth",
    verbose: bool = False,
) -> dict:
    """_summary_

    Args:
        bam (str): Path to BAM file
        bai (str, optional): Path to BAM index file. Defaults to None.
        mt_ref (str, optional): Mitochondrial reference. Defaults to "rcrs".
        reference_fa (str, optional): Referenoe genome. Defaults to None.
        contig_name (str, optional): Name of mitochondrial contig. Defaults to None.
        out_dir (str, optional): Output directory. Defaults to None.
        tmp_dir (str, optional): Tmp directory. Defaults to None.
        remove_tmp (bool, optional): Remove tmp directory. Defaults to False.
        prefix (str, optional): Prefix. Defaults to None.
        ncores (int, optional): Number of cores. Defaults to 1.
        m2_extra_args (str, optional): Extra args for Mutect2. Defaults to None.
        f_score_beta (float, optional): F score beta. Defaults to 1.
        contamination_filter (bool, optional): Contamination filter. Defaults to True.
        max_alt_allele_count (int, optional): Maximuam alt allele count. Defaults to 4.
        vaf_treshold (float, optional): Minimum variant allele fraction. Defaults to 0.
        autosomal_coverage (float, optional): median autosomal coverge. Defaults to 0.
        blacklisted_sites (str, optional): Custom BED file with blacklisted sites. Defaults to None.
        remove_non_pass (bool, optional): Remove nonpass variants. Defaults to True.
        normalize (bool, optional): Normalize variants. Defaults to True.
        min_hom_treshold (float, optional): Minimal homoplasmy treshold. Defaults to 0.95.
        split_strands (bool, optional): Split strands of mitochondrial genome in visualziation. Defaults to True.
        population_freqs (bool, optional): Add population frequencies. Defaults to True.
        patho_predictions (bool, optional): Add pathogenicity predictions. Defaults to True.
        phenotype_annot (bool, optional): Add phenotype annotations. Defaults to True.
        conservation_scores (bool, optional): Add conservation scores. Defaults to True.
        save_as_png (bool, optional): Save vis plot as PNG. Defaults to True.
        create_annotation_report (bool, optional): Create CSV annotation report. Defaults to True.
        snpeff_path (str, optional): Path to SnpEff. Defaults to "snpeff".
        snpsift_path (str, optional): Path to SnpSift. Defaults to "snpsift".
        haplogrep3_path (str, optional): Path to Haplogrep3. Defaults to "haplogrep3".
        gatk_path (str, optional): Path to GATK. Defaults to "gatk".
        bwamem2_path (str, optional): Path to bwa-mem2. Defaults to "bwa-mem2".
        haplocheck_path (str, optional): Path to Haplocheck. Defaults to "haplocheck".
        mosdepth_path (str, optional): Path to mosdepth. Defaults to "mosdepth".
        verbose (bool, optional): Verbosity. Defaults to False.

    Returns:
        dict: Main output file paths
    """
    if not prefix:
        prefix = get_file_basename(bam)

    if not out_dir:
        out_dir = get_file_directory(bam)

    if not tmp_dir:
        tmp_dir = "tmp"

    results_dir = "results"
    intermediates = os.path.join(out_dir, tmp_dir)
    results = os.path.join(out_dir, results_dir)

    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(intermediates, exist_ok=True)
    os.makedirs(results, exist_ok=True)

    final_outputs = {}

    logging.info(f"Running mitopy pipeline on {bam} file.")

    # Preprocess BAM file
    prep_out = do_preprocess(
        bam=bam,
        bai=bai,
        reference_fa=reference_fa,
        contig_name=contig_name,
        out_dir=f"{intermediates}/preprocess",
        prefix=prefix,
        gatk_path=gatk_path,
        verbose=verbose,
    )

    # Align to canonical reference
    align_canonical = do_align(
        ubam=prep_out["unmapped_bam"],
        mt_ref=mt_ref,
        out_dir=f"{intermediates}/align",
        prefix=prefix,
        gatk_path=gatk_path,
        bwamem2_path=bwamem2_path,
        ncores=ncores,
        verbose=verbose,
    )

    final_outputs.update(align_canonical)

    # Align to shifted reference
    align_shifted = do_align(
        ubam=prep_out["unmapped_bam"],
        mt_ref=mt_ref,
        out_dir=f"{intermediates}/align",
        prefix=prefix,
        gatk_path=gatk_path,
        bwamem2_path=bwamem2_path,
        ncores=ncores,
        shifted=True,
        verbose=verbose,
    )

    final_outputs.update(
        {"shifted_" + key: value for key, value in align_shifted.items()}
    )

    # Call variants in non-control region
    call_canonical = do_call(
        bam=align_canonical["dedup_sorted_bam"],
        mt_ref=mt_ref,
        m2_extra_args=m2_extra_args,
        out_dir=f"{intermediates}/call",
        prefix=prefix,
        gatk_path=gatk_path,
        verbose=verbose,
    )

    # Call variants in control region
    call_shifted = do_call(
        bam=align_shifted["dedup_sorted_bam"],
        mt_ref=mt_ref,
        m2_extra_args=m2_extra_args,
        out_dir=f"{intermediates}/call",
        prefix=prefix,
        gatk_path=gatk_path,
        shifted=True,
        verbose=verbose,
    )

    # Merge variant calls
    merged = do_merge(
        vcf=call_canonical["raw_vcf"],
        vcf_shifted=call_shifted["raw_vcf"],
        mt_ref=mt_ref,
        stats=call_canonical["raw_vcf_stats"],
        stats_shifted=call_shifted["raw_vcf_stats"],
        out_dir=f"{intermediates}/merge",
        prefix=prefix,
        gatk_path=gatk_path,
        verbose=verbose,
    )

    # Postprocessing
    postprocessed = do_postprocess(
        vcf=merged["merged_vcf"],
        stats=merged["merged_vcf_stats"],
        mt_ref=mt_ref,
        f_score_beta=f_score_beta,
        contamination_filter=contamination_filter,
        max_alt_allele_count=max_alt_allele_count,
        vaf_treshold=vaf_treshold,
        autosomal_coverage=autosomal_coverage,
        blacklisted_sites=blacklisted_sites,
        remove_non_pass=remove_non_pass,
        normalize=normalize,
        out_dir=f"{intermediates}/postprocess",
        prefix=prefix,
        gatk_path=gatk_path,
        haplocheck_path=haplocheck_path,
        verbose=verbose,
    )

    final_outputs.update(postprocessed)

    # Calculate coverage
    coverage = do_coverage(
        mt_bam=align_canonical["dedup_sorted_bam"],
        shifted_mt_bam=align_shifted["dedup_sorted_bam"],
        prefix=prefix,
        out_dir=f"{intermediates}/coverage",
        mosdepth_path=mosdepth_path,
        verbose=verbose,
    )

    # Annotate
    annotation = do_annotate(
        vcf=postprocessed["postprocessed_vcf"],
        min_hom_treshold=min_hom_treshold,
        population_freqs=population_freqs,
        patho_predictions=patho_predictions,
        phenotype_annot=phenotype_annot,
        conservation_scores=conservation_scores,
        prefix=prefix,
        create_csv=create_annotation_report,
        out_dir=f"{intermediates}/annotate",
        snpeff_path=snpeff_path,
        snpsift_path=snpsift_path,
        verbose=verbose,
    )

    final_outputs.update(annotation)

    # Visualize
    visualization = do_visualize(
        vcf=postprocessed["postprocessed_vcf"],
        coverage_csv=coverage["coverage_csv"],
        split_strands=split_strands,
        save_as_png=save_as_png,
        out_dir=f"{intermediates}/visualize",
        prefix=prefix,
    )

    final_outputs.update(visualization)

    # Identify haplogroups
    haplogroups = do_identify_haplogroup(
        vcf=postprocessed["postprocessed_vcf"],
        mt_ref=mt_ref,
        haplogrep3_path=haplogrep3_path,
        prefix=prefix,
        out_dir=f"{intermediates}/haplogroup",
        verbose=verbose,
    )

    final_outputs.update(haplogroups)

    # Move final outputs to results directory
    for _, fpath in final_outputs.items():
        shutil.copy(fpath, results)

    # Remove intermediates
    if remove_tmp:
        shutil.rmtree(intermediates)

    # Check output files
    if check_files_exist(list(final_outputs.values())):
        logging.info(
            f"Mitopy pipeline completed successfully! See the final results: {results}"
        )
    else:
        logging.error("Some output files are missing! Please rerun the analysis.")
        sys.exit(1)

    return final_outputs
