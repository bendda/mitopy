import os
import sys

MITOPY_DIR = os.path.dirname(sys.modules["mitopy"].__file__)
DATA_DIR = os.path.join(MITOPY_DIR, "data")

REF_DIR = os.path.join(DATA_DIR, "mt_reference")
FILTER_DIR = os.path.join(DATA_DIR, "filter_data")
ANNOT_DIR = os.path.join(DATA_DIR, "annotation_data")
VIS_DIR = os.path.join(DATA_DIR, "vis_data")


MT_REFS = {
    "rcrs": f"{REF_DIR}/rcrs/rcrs.fasta",
    "rcrs_shifted": f"{REF_DIR}/rcrs/rcrs_shifted.fasta",
    "rcrs_shift_back_chain": f"{REF_DIR}/rcrs/rcrs.shift_back.chain",
    "rsrs": f"{REF_DIR}/rsrs/rsrs.fasta",
    "rsrs_shifted": f"{REF_DIR}/rsrs/rsrs_shifted.fasta",
    "rsrs_shift_back_chain": f"{REF_DIR}/rsrs/rsrs.shift_back.chain",
}


MT_BLACKLIST = {
    "rcrs": f"{FILTER_DIR}/rcrs_blacklist.bed",
    "rsrs": f"{FILTER_DIR}/rsrs_blacklist.bed",
}


ANNOT_RESOURCES = {
    "mitotip": f"{ANNOT_DIR}/mitotip.csv",
    "pon_mt_trna": f"{ANNOT_DIR}/pon_mt_trna.csv",
    "gnomad": f"{ANNOT_DIR}/gnomad_freq.csv",
    "mitomap": f"{ANNOT_DIR}/mitomap.csv",
    "general": f"{ANNOT_DIR}/general_annot.csv",
    "snpeff_config": f"{ANNOT_DIR}/snpeff/snpeff.config",
    "snpeff_db": "NC_012920",
    "conservation_scores": f"{ANNOT_DIR}/conservation_scores.csv",
    "clinvar": f"{ANNOT_DIR}/clinvar.csv",
    "dbsnp": f"{ANNOT_DIR}/dbsnp.csv",
    "sift": f"{ANNOT_DIR}/sift.csv",
}


VIS_RESOURCES = {"mitomap": f"{VIS_DIR}/mitomap.csv"}
