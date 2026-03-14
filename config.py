"""
config.py — Central configuration for the COVID-19 Mutation Tracker
All paths, reference sequences, variant definitions, and colour palettes live here.
"""

import os
from pathlib import Path

# ── Paths ──────────────────────────────────────────────────────────────────
BASE_DIR   = Path(__file__).parent
DATA_DIR   = BASE_DIR / "data"
OUTPUT_DIR = BASE_DIR / "data" / "processed"
RAW_DIR    = BASE_DIR / "data" / "raw"

for d in [DATA_DIR, OUTPUT_DIR, RAW_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# ── Reference genome ───────────────────────────────────────────────────────
REFERENCE_ACCESSION = "MN908947"          # Wuhan-Hu-1 (NCBI)
REFERENCE_FASTA     = RAW_DIR / "MN908947_reference.fasta"
REFERENCE_GFF       = RAW_DIR / "MN908947_annotation.gff"

# ── Spike protein (gene S) coordinates on MN908947 ────────────────────────
SPIKE_START = 21562   # 1-based, inclusive
SPIKE_END   = 25384
SPIKE_AA_LEN = 1273

# ── Spike domain boundaries (amino-acid positions, 1-based) ───────────────
SPIKE_DOMAINS = {
    "Signal peptide": (1,   13),
    "NTD":            (14,  305),
    "RBD":            (306, 527),
    "RBM":            (438, 506),
    "Central domain": (528, 680),
    "Furin cleavage": (675, 695),
    "Fusion peptide": (788, 806),
    "HR1":            (912, 984),
    "HR2":            (1163,1213),
    "TM":             (1213,1237),
    "Cytoplasmic":    (1237,1273),
}

# ── Known variant lineage definitions ─────────────────────────────────────
# Defining spike mutations for each WHO-labelled variant
VARIANT_DEFINITIONS = {
    "Alpha (B.1.1.7)": {
        "pango": "B.1.1.7",
        "who":   "Alpha",
        "defining_spike": ["N501Y", "P681H", "D614G", "HV69-70del", "Y144del"],
        "origin_country": "United Kingdom",
        "first_detected": "2020-09",
        "colour": "#1D9E75",
    },
    "Beta (B.1.351)": {
        "pango": "B.1.351",
        "who":   "Beta",
        "defining_spike": ["K417N", "E484K", "N501Y", "D614G"],
        "origin_country": "South Africa",
        "first_detected": "2020-10",
        "colour": "#3B8BD4",
    },
    "Gamma (P.1)": {
        "pango": "P.1",
        "who":   "Gamma",
        "defining_spike": ["K417T", "E484K", "N501Y", "D614G"],
        "origin_country": "Brazil",
        "first_detected": "2020-11",
        "colour": "#7F77DD",
    },
    "Delta (B.1.617.2)": {
        "pango": "B.1.617.2",
        "who":   "Delta",
        "defining_spike": ["L452R", "T478K", "P681R", "D614G"],
        "origin_country": "India",
        "first_detected": "2021-04",
        "colour": "#EF9F27",
    },
    "Omicron BA.1": {
        "pango": "B.1.1.529",
        "who":   "Omicron",
        "defining_spike": ["K417N", "E484A", "N501Y", "D614G", "P681H",
                           "G339D", "S371L", "S373P", "S375F", "Q493R",
                           "G496S", "Q498R", "Y505H", "T547K", "H655Y",
                           "N679K", "N764K", "D796Y", "N856K", "Q954H",
                           "N969K", "L981F"],
        "origin_country": "South Africa",
        "first_detected": "2021-11",
        "colour": "#D4537E",
    },
    "Omicron BA.2": {
        "pango": "BA.2",
        "who":   "Omicron",
        "defining_spike": ["K417N", "N501Y", "D614G", "T478K", "G339D",
                           "S371F", "S373P", "S375F", "T376A", "D405N",
                           "R408S", "Q493R", "Q498R", "Y505H"],
        "origin_country": "Multiple",
        "first_detected": "2021-12",
        "colour": "#993556",
    },
    "XBB.1.5": {
        "pango": "XBB.1.5",
        "who":   "Omicron (recombinant)",
        "defining_spike": ["G339H", "R346T", "L368I", "V445P", "G446S",
                           "N460K", "F486P", "F490S", "R493Q"],
        "origin_country": "United States",
        "first_detected": "2022-10",
        "colour": "#A32D2D",
    },
    "JN.1": {
        "pango": "JN.1",
        "who":   "Omicron (BA.2.86 sub)",
        "defining_spike": ["L455S", "F456L", "K417N", "N501Y", "D614G",
                           "P681H", "N460K", "R346T"],
        "origin_country": "Denmark",
        "first_detected": "2023-08",
        "colour": "#71243E",
    },
}

# ── Key mutations to track ─────────────────────────────────────────────────
KEY_MUTATIONS = [
    # mutation, protein position, functional category, notes
    ("D614G",  614, "Transmissibility",  "Universal backbone; increased virion stability"),
    ("N501Y",  501, "ACE2 affinity",     "Enhanced ACE2 binding; present Alpha/Beta/Gamma/Omicron"),
    ("E484K",  484, "Immune escape",     "Antibody neutralisation evasion; convergent evolution"),
    ("E484A",  484, "Immune escape",     "Omicron-specific E484 substitution"),
    ("K417N",  417, "Immune escape",     "Reduces antibody binding; Beta/Omicron"),
    ("K417T",  417, "Immune escape",     "Gamma-specific variant of K417 position"),
    ("L452R",  452, "Infectivity",       "Increased ACE2 affinity; Delta dominant"),
    ("T478K",  478, "ACE2 binding",      "Charge change at ACE2 interface; Delta/Omicron"),
    ("P681H",  681, "Furin cleavage",    "Furin site enhancement; Alpha/Omicron"),
    ("P681R",  681, "Furin cleavage",    "Stronger furin cleavage; Delta"),
    ("H655Y",  655, "Furin cleavage",    "Upstream furin site; Gamma/Omicron"),
    ("N679K",  679, "Furin cleavage",    "Additional positive charge; Omicron BA.1"),
    ("G339D",  339, "RBD stability",     "Structural; Omicron clade"),
    ("R346T",  346, "Immune escape",     "Monoclonal antibody escape; XBB/JN.1"),
    ("L455S",  455, "ACE2 binding",      "JN.1 defining; ACE2 affinity change"),
    ("F486P",  486, "ACE2/antibody",     "XBB.1.5 key substitution"),
]

# ── Colour palette ─────────────────────────────────────────────────────────
FUNCTIONAL_COLOURS = {
    "Transmissibility": "#3B8BD4",
    "ACE2 affinity":    "#1D9E75",
    "Immune escape":    "#D4537E",
    "Furin cleavage":   "#EF9F27",
    "RBD stability":    "#7F77DD",
    "ACE2 binding":     "#5DCAA5",
    "ACE2/antibody":    "#A32D2D",
}

VARIANT_COLOURS = {v["who"]: v["colour"]
                   for v in VARIANT_DEFINITIONS.values()}

# ── Nextstrain metadata URL (public, no login required) ────────────────────
NEXTSTRAIN_METADATA_URL = (
    "https://data.nextstrain.org/files/ncov/open/metadata.tsv.zst"
)
NEXTSTRAIN_SEQUENCES_URL = (
    "https://data.nextstrain.org/files/ncov/open/sequences.fasta.zst"
)
