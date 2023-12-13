import logging
import pandas as pd
import vcfpy
import numpy as np
from .utils import (
    get_file_basename,
    create_output_path,
    get_file_directory,
    check_files_exist,
)
import os
import sys
from .executable import Executable
from .constants import ANNOT_RESOURCES


general_fields = {
    "LOCUS": [
        ("ID", "LOCUS"),
        ("Number", "1"),
        ("Type", "String"),
        ("Description", "Variant locus/gene"),
    ],
    "BIOTYPE": [
        ("ID", "BIOTYPE"),
        ("Number", "1"),
        ("Type", "String"),
        ("Description", "Locus/gene biotype"),
    ],
}

conservation_fields = {
    "phastCons100way": [
        ("ID", "phastCons100way"),
        ("Number", "1"),
        ("Type", "Float"),
        (
            "Description",
            "PhastCons conservation score (conserved: score > 0.7 [soft treshold])",
        ),
    ],
    "phyloP100way": [
        ("ID", "phyloP100way"),
        ("Number", "1"),
        ("Type", "Float"),
        ("Description", "PhyloP conservation score (conserved: score > 0)"),
    ],
}


clinvar_fields = {
    "ClinVar_ID": [
        ("ID", "ClinVar_ID"),
        ("Number", "1"),
        ("Type", "String"),
        (
            "Description",
            "ClinVar Variation ID",
        ),
    ],
    "CLNDN": [
        ("ID", "CLNDN"),
        ("Number", "1"),
        ("Type", "String"),
        ("Description", "ClinVar disease name associated with variant"),
    ],
    "CLNSIG": [
        ("ID", "CLNSIG"),
        ("Number", "1"),
        ("Type", "String"),
        ("Description", "Clinical significance for variant"),
    ],
    "CLNDISDB": [
        ("ID", "CLNDISDB"),
        ("Number", "1"),
        ("Type", "String"),
        ("Description", "Disease database name and identifier for variant"),
    ],
}


gnomad_fields = {
    "GNOMAD_AC_HOM": [
        ("ID", "GNOMAD_AC_HOM"),
        ("Number", "1"),
        ("Type", "Float"),
        (
            "Description",
            "Gnomad allele count restricted to variants with a heteroplasmy level >= 0.95",
        ),
    ],
    "GNOMAD_AF_HOM": [
        ("ID", "GNOMAD_AF_HOM"),
        ("Number", "1"),
        ("Type", "Float"),
        (
            "Description",
            "Gnomad allele frequency restricted to variants with a heteroplasmy level >= 0.95",
        ),
    ],
    "GNOMAD_AF_HET": [
        ("ID", "GNOMAD_AF_HET"),
        ("Number", "1"),
        ("Type", "Float"),
        (
            "Description",
            "Gnomad allele frequency restricted to variants with a heteroplasmy level >= 0.10 and < 0.95",
        ),
    ],
    "GNOMAD_AC_HET": [
        ("ID", "GNOMAD_AC_HET"),
        ("Number", "1"),
        ("Type", "Float"),
        (
            "Description",
            "Gnomad allele count restricted to variants with a heteroplasmy level >= 0.10 and < 0.95",
        ),
    ],
}


mitotip_fields = {
    "MitoTIP_Score": [
        ("ID", "MitoTIP_Score"),
        ("Number", "1"),
        ("Type", "Float"),
        ("Description", "tRNA raw pathogenicity score from Mitotip"),
    ],
    "MitoTIP_Prediction": [
        ("ID", "MitoTIP_Prediction"),
        ("Number", "1"),
        ("Type", "String"),
        ("Description", "tRNA pathogenicity classification from Mitotip"),
    ],
}


pontrna_fields = {
    "PONmttRNA_Probability": [
        ("ID", "PONmttRNA_Probability"),
        ("Number", "1"),
        ("Type", "String"),
        ("Description", "tRNA probability of pathogenicity from PON-mt-tRNA"),
    ],
    "PONmttRNA_Prediction": [
        ("ID", "PONmttRNA_Prediction"),
        ("Number", "1"),
        ("Type", "Float"),
        ("Description", "tRNA pathogenicity classification from PON-mt-tRNA"),
    ],
}


sift_fields = {
    "SIFT": [
        ("ID", "SIFT"),
        ("Number", "1"),
        ("Type", "String"),
        ("Description", "Pathogenicity classification from SIFT"),
    ],
    "SIFT_score": [
        ("ID", "SIFT_score"),
        ("Number", "1"),
        ("Type", "Float"),
        ("Description", "Pathogenicity prediction score from SIFT"),
    ],
}


mitomap_fields = {
    "MITOMAP_GENBANK_AC": [
        ("ID", "MITOMAP_GENBANK_AC"),
        ("Number", "1"),
        ("Type", "Float"),
        (
            "Description",
            "Allele count in GenBank out of 61168 full length (FL) human chrM sequences",
        ),
    ],
    "MITOMAP_GENBANK_AF": [
        ("ID", "MITOMAP_GENBANK_AF"),
        ("Number", "1"),
        ("Type", "Float"),
        ("Description", "Allele Frequency in full length (FL) Genbank sequence set"),
    ],
    "MITOMAP_PubmedIDs": [
        ("ID", "MITOMAP_PubmedIDs"),
        ("Number", "1"),
        ("Type", "Integer"),
        ("Description", "Pubmed IDs"),
    ],
    "MITOMAP_Disease": [
        ("ID", "MITOMAP_Disease"),
        ("Number", "1"),
        ("Type", "String"),
        ("Description", "Putative Disease Association"),
    ],
    "MITOMAP_DiseaseStatus": [
        ("ID", "MITOMAP_DiseaseStatus"),
        ("Number", "1"),
        ("Type", "String"),
        ("Description", "Disease Association Status"),
    ],
}


def _snpeff_annotate(in_vcf: str, out_vcf: str, snpeff_exec: Executable) -> None:
    """Perform functional annotation using SnpEff."""

    db = ANNOT_RESOURCES["snpeff_db"]

    params = {
        "-c": ANNOT_RESOURCES["snpeff_config"],
        "-noStats": True,
        "-no-upstream": True,
        "-no-downstream": True,
    }

    snpeff_exec.run(db, in_vcf, redirect_out=out_vcf, **params)


def do_annotate(
    vcf: str,
    min_hom_treshold: float = 0.95,
    population_freqs: bool = True,
    patho_predictions: bool = True,
    phenotype_annot: bool = True,
    conservation_scores: bool = True,
    prefix: str = None,
    out_dir: str = None,
    create_csv: bool = True,
    snpeff_path: str = "snpeff",
    snpsift_path: str = "snpsift",
    verbose: bool = False,
) -> dict:
    """Annotate mitochondrial variants.

    Args:
        vcf (str): Path to input VCF file
        min_hom_treshold (float, optional): Minimum homoplasmy level treshold. Defaults to 0.95.
        population_freqs (bool, optional): Annotate variants with population frequencies. Defaults to True.
        patho_predictions (bool, optional): Annotate variants with in-silico pathogenicity predictions. Defaults to True.
        phenotype_annot (bool, optional): Annotate variants with phenotype information. Defaults to True.
        conservation_scores (bool, optional): Annotate with conservation scores. Defaults to True.
        prefix (str, optional): Prefix. Defaults to None.
        out_dir (str, optional): Output directory. Defaults to None.
        create_csv (bool, optional): Export annotated variants to CSV format. Defaults to True.
        snpeff_path (str, optional): Path to snpeff executable. Defaults to "snpeff".
        snpsift_path (str, optional): Path to snpsift executable. Defaults to "snpsift".
        verbose (bool, optional): Verbosity. Defaults to False.

    Returns:
        dict: Main output file paths
    """
    snpeff = Executable(snpeff_path, verbose)

    if not prefix:
        prefix = get_file_basename(vcf)

    if not out_dir:
        out_dir = get_file_directory(vcf)

    os.makedirs(out_dir, exist_ok=True)

    out_vcf = create_output_path(prefix, out_dir, "_annotated", ".vcf")
    snpeff_vcf = create_output_path(prefix, out_dir, "_snpeff", ".vcf")

    # Perform functional annotation
    logging.info("Performing functional annotation using SNPeff...")
    _snpeff_annotate(in_vcf=vcf, out_vcf=snpeff_vcf, snpeff_exec=snpeff)

    # Variant specific annotation db
    info_fields = {}
    info_df = pd.DataFrame(columns=["POS", "REF", "ALT"])

    # Site specific annotation db
    site_fields = {}
    site_annot_df = pd.DataFrame(columns=["POS"])

    # Load general annotation
    site_fields.update(general_fields)
    general_df = pd.read_csv(ANNOT_RESOURCES["general"])
    site_annot_df = site_annot_df.merge(general_df, how="outer", on=["POS"])

    # Load conservation scores
    if conservation_scores:
        logging.info("Adding conservation scores...")
        site_fields.update(conservation_fields)
        cons_df = pd.read_csv(ANNOT_RESOURCES["conservation_scores"])
        site_annot_df = site_annot_df.merge(cons_df, how="outer", on=["POS"])

    # Load pathogenicity predictions
    if patho_predictions:
        logging.info("Adding pathogenicity predictions...")
        # Add SIFT
        info_fields.update(sift_fields)
        patho_df = pd.read_csv(ANNOT_RESOURCES["sift"])
        info_df = info_df.merge(patho_df, how="outer", on=["POS", "REF", "ALT"])

        # Add mitotip
        info_fields.update(mitotip_fields)
        mitotip_df = pd.read_csv(ANNOT_RESOURCES["mitotip"])
        info_df = info_df.merge(mitotip_df, how="outer", on=["POS", "REF", "ALT"])

        # Add PON-mt-trna
        info_fields.update(pontrna_fields)
        pontrna_df = pd.read_csv(ANNOT_RESOURCES["pon_mt_trna"])
        info_df = info_df.merge(pontrna_df, how="outer", on=["POS", "REF", "ALT"])

    # Load population frequencies
    if population_freqs:
        logging.info("Adding population frequencies...")

        # Add gnomad
        info_fields.update(gnomad_fields)
        gnomad_df = pd.read_csv(ANNOT_RESOURCES["gnomad"])
        info_df = info_df.merge(gnomad_df, how="outer", on=["POS", "REF", "ALT"])

    # Load phenotype annotations
    if phenotype_annot:
        logging.info("Adding phenotype annotations...")

        # Add MITOMAP
        info_fields.update(mitomap_fields)
        mitomap_df = pd.read_csv(ANNOT_RESOURCES["mitomap"])
        info_df = info_df.merge(mitomap_df, how="outer", on=["POS", "REF", "ALT"])

        # Add ClinVat
        info_fields.update(clinvar_fields)
        clinvar_df = pd.read_csv(ANNOT_RESOURCES["clinvar"])
        info_df = info_df.merge(clinvar_df, how="outer", on=["POS", "REF", "ALT"])

    # Create search database
    info_df = info_df.replace({np.nan: None})
    info_df["key"] = info_df.apply(
        lambda row: (row["POS"], row["REF"], row["ALT"]), axis=1
    )

    reader = vcfpy.Reader.from_path(snpeff_vcf)

    # Add headers
    for header in {**info_fields, **site_fields}.values():
        reader.header.add_info_line(vcfpy.OrderedDict(header))

    logging.info("Writing additional annotations...")
    # Write output vcf
    writer = vcfpy.Writer.from_path(out_vcf, reader.header)

    for record in reader:
        # Add site annotations
        site_annot = site_annot_df.query("POS == @record.POS")
        if not site_annot.empty:
            for field in site_fields:
                record.INFO[field] = site_annot[field].values[0]

        # Add variant annotations
        rec_key = (record.POS, record.REF, record.ALT[0].serialize())
        annot = info_df.query("key == @rec_key")
        if not annot.empty:
            for field in info_fields:
                if annot[field].values[0]:
                    record.INFO[field] = annot[field].values[0]

        # Adjust genotype based on min_hom_treshold
        gt_info = record.calls[0]
        genotype = "1/1" if gt_info.data["AF"][0] > min_hom_treshold else "0/1"
        gt_info.set_genotype(genotype)

        # Write the updated record to the output VCF file
        writer.write_record(record)

    reader.close()
    writer.close()

    output_paths = {
        "annotated_vcf": out_vcf,
    }

    # Create CSV report
    if create_csv:
        logging.info("Creating CSV report...")
        csv_fields = []

        variant_info_fields = list(info_fields.keys())
        site_info_fields = list(site_fields.keys())

        standard_fields = ["CHROM", "POS", "REF", "ALT"]
        genotype_fields = ["GEN[*].AF", "GEN[*].GT"]
        snpeff_fields = [
            "ANN[*].EFFECT",
            "ANN[*].IMPACT",
            "ANN[*].GENE",
            "ANN[*].GENEID",
            "ANN[*].FEATURE",
            "ANN[*].FEATUREID",
            "ANN[*].BIOTYPE",
            "ANN[*].HGVS_C",
            "ANN[*].HGVS_P",
        ]

        csv_fields.extend(standard_fields)
        csv_fields.extend(site_info_fields)
        csv_fields.extend(genotype_fields)
        csv_fields.extend(snpeff_fields)
        csv_fields.extend(variant_info_fields)

        out_csv = create_output_path(prefix, out_dir, "_annotated", ".csv")

        logging.info("Extracting relevant fields from VCF...")
        snpsift = Executable(snpsift_path, verbose)

        snpsift.run(
            out_vcf,
            " ".join(csv_fields),
            subcommand="extractFields",
            redirect_out=out_csv,
            **{"-e": "."},
        )

        # Postprocess CSV
        out_df = pd.read_csv(out_csv, sep="\t")
        out_df.columns = out_df.columns.map(
            lambda x: "snpEff_" + x.split(".")[1] if x.startswith("ANN[*]") else x
        )
        out_df.rename(columns={"GEN[*].AF": "Heteroplasmy Fraction"}, inplace=True)
        out_df.rename(columns={"GEN[*].GT": "MT Variant Type"}, inplace=True)

        # Annotate variant types
        out_df["MT Variant Type"] = [
            "homoplasmic" if gt == "1/1" else "heteroplasmic"
            for gt in out_df["MT Variant Type"]
        ]
        out_df.to_csv(out_csv, index=False)
        output_paths["annot_csv"] = out_csv

    # Check if output files exist
    if check_files_exist(list(output_paths.values())):
        logging.info(f"Annotation completed successfully.")
    else:
        logging.error("Some output files are missing! Please rerun the analysis.")
        sys.exit(1)

    return output_paths
