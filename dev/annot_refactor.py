import pandas as pd
import subprocess

# Prepare annotation files


annot_resources = {
    "mitotip": "mitotip_scores.txt",
    "pon_mt_trna": "PON-mt-tRNA_predictions.txt",
    "gnomad": "gnomad.genomes.v3.1.sites.chrM.reduced_annotations.tsv",
    "mitomap": "mitomap_disease.vcf",
    "clinvar": "clinvar.vcf.gz",
    "phast": "chrM.phastCons100way.wigFix",
    "phylo": "chrM.phyloP100way.wigFix",
    "sift": "MT.gz",
}


def refactor_mitotip():
    """Refactor MITOTIP annotaion file."""
    mitotip_df = pd.read_csv(annot_resources["mitotip"], sep="\t")
    mitotip_df.rename(
        columns={"Position": "POS", "rCRS": "REF", "Alt": "ALT"}, inplace=True
    )

    # Score interpretation
    quartile_interpretation = {
        "Q1": "likely pathogenic",
        "Q2": "possibly pathogenic",
        "Q3": "possibly benign",
        "Q4": "likely benign",
    }

    mitotip_df["MitoTIP_Prediction"] = [
        quartile_interpretation[q] for q in mitotip_df["Quartile"]
    ]
    mitotip_df.drop(
        columns=["Quartile", "Count", "Percentage", "Mitomap_Status"], inplace=True
    )
    mitotip_df.to_csv("mitotip.csv", index=False)


def refactor_pontrna():
    """Refactor PON-mt-trna annotation file."""
    pon_mt_trna_df = pd.read_csv(annot_resources["pon_mt_trna"], sep="\t", skiprows=25)
    pon_mt_trna_df.rename(
        columns={
            "mtDNA_position": "POS",
            "Reference_nucleotide": "REF",
            "New_nucleotide": "ALT",
            "Classification": "PONmttRNA_Prediction",
            "ML_probability_of_pathogenicity": "PONmttRNA_Probability",
        },
        inplace=True,
    )
    pon_mt_trna_df.drop(
        columns=["mt-tRNA", "Prediction_of_Kondrashov_method"], inplace=True
    )
    pon_mt_trna_df.to_csv("pon_mt_trna.csv", index=False)


def refactor_gnomad():
    """Refactor gnomad annotation file."""
    gnomad = pd.read_csv(annot_resources["gnomad"], sep="\t")
    gnomad.rename(
        columns={
            "position": "POS",
            "ref": "REF",
            "alt": "ALT",
            "AC_het": "GNOMAD_AC_HET",
            "AC_hom": "GNOMAD_AC_HOM",
            "AF_het": "GNOMAD_AF_HET",
            "AF_hom": "GNOMAD_AF_HOM",
        },
        inplace=True,
    )
    gnomad.drop(
        columns=["chromosome", "filters", "AN", "max_observed_heteroplasmy"],
        inplace=True,
    )
    gnomad.to_csv("gnomad_freq.csv", index=False)


def refactor_mitomap():
    """Refactor mitomap annotation file."""

    # normalize and convert to TSV
    subprocess.run(
        f"bcftools norm -m - -o mitomap_norm.vcf {annot_resources['mitomap']}",
        shell=True,
    )
    subprocess.run(
        f"gatk VariantsToTable -O mitomap.tsv -V {annot_resources['mitomap']}"
    )

    mitomap = pd.read_csv("mitomap.tsv", sep="\t")
    mitomap.rename(columns={"AC": "GENBANK_AC", "AF": "GENBANK_AF"}, inplace=True)
    mitomap.drop(columns=["ID", "CHROM", "QUAL", "FILTER", "HGFL"], inplace=True)
    mitomap.columns = [
        "MITOMAP_" + col if col not in ["POS", "REF", "ALT"] else col
        for col in mitomap.columns
    ]
    mitomap.to_csv("mitomap.csv", index=False)


def refactor_clinvar():
    """Refactor Clinvar annotation file."""
    subprocess.run(
        f"bcftools view {annot_resources['clinvar']} --regions MT > clinvar_MT.tsv",
        shell=True,
    )
    clinvar = pd.read_csv("clinvar_MT.tsv", sep="\t")
    clinvar = clinvar[["POS", "ID", "REF", "ALT", "CLNDN", "CLNSIG", "CLNDISDB"]]
    clinvar.rename(columns={"ID": "ClinVar_ID"}, inplace=True)
    clinvar.to_csv("clinvar.csv", index=False)


def refactor_conservation_scores():
    """Refactor conservation scores."""
    phast = pd.read_csv(annot_resources["phast"], names=["phastCons100way"], skiprows=1)
    phylo = pd.read_csv(annot_resources["phylo"], names=["phyloP100way"], skiprows=1)

    cons = pd.DataFrame()
    cons["POS"] = range(1, 16570)

    cons = pd.concat([cons, phast, phylo], axis=1)
    cons.to_csv("conservation_scores.csv", index=False)


def refactor_sift():
    """Refactor sift annotation file."""
    sift = pd.read_csv(annot_resources["sift"], sep="\t")
    sift.rename(
        columns={"Position": "POS", "Ref_allele": "REF", "New_allele": "ALT"},
        inplace=True,
    )
    sift["SIFT"] = [
        "neutral" if score > 0.05 else "deleterious" for score in sift["SIFT_score"]
    ]

    sift = sift[["POS", "REF", "ALT", "SIFT", "SIFT_score"]]
    sift.to_csv("sift.csv", index=False)
