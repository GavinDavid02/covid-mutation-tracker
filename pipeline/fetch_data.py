"""
pipeline/fetch_data.py
──────────────────────
Downloads:
  1. Wuhan-Hu-1 reference genome (MN908947) from NCBI Entrez
  2. Nextstrain open metadata (TSV, zstd-compressed, ~500 MB)
  3. Optionally a GISAID FASTA if you supply a local path

Run:
    python pipeline/fetch_data.py
    python pipeline/fetch_data.py --source nextstrain --max-sequences 5000
    python pipeline/fetch_data.py --source gisaid --gisaid-fasta /path/to/sequences.fasta
"""

import argparse
import gzip
import io
import os
import shutil
import sys
import time
import urllib.request
from pathlib import Path

import pandas as pd
from Bio import Entrez, SeqIO
from tqdm import tqdm

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))
from config import (
    RAW_DIR, OUTPUT_DIR,
    REFERENCE_ACCESSION, REFERENCE_FASTA,
    NEXTSTRAIN_METADATA_URL,
)

Entrez.email = "researcher@example.com"   # NCBI requires an e-mail


# ── Helpers ────────────────────────────────────────────────────────────────

class DownloadProgressBar(tqdm):
    def update_to(self, b=1, bsize=1, tsize=None):
        if tsize is not None:
            self.total = tsize
        self.update(b * bsize - self.n)


def download_url(url: str, dest: Path, desc: str = ""):
    """Download *url* to *dest* with a tqdm progress bar."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    with DownloadProgressBar(unit="B", unit_scale=True, miniters=1, desc=desc) as t:
        urllib.request.urlretrieve(url, filename=dest, reporthook=t.update_to)


# ── Reference genome ───────────────────────────────────────────────────────

def fetch_reference_genome() -> Path:
    """Download MN908947 (Wuhan-Hu-1) from NCBI if not already present."""
    if REFERENCE_FASTA.exists():
        print(f"[✓] Reference already present: {REFERENCE_FASTA}")
        return REFERENCE_FASTA

    print(f"[→] Fetching reference {REFERENCE_ACCESSION} from NCBI …")
    retries = 3
    for attempt in range(retries):
        try:
            handle = Entrez.efetch(
                db="nucleotide",
                id=REFERENCE_ACCESSION,
                rettype="fasta",
                retmode="text",
            )
            record = SeqIO.read(handle, "fasta")
            handle.close()
            with open(REFERENCE_FASTA, "w") as fh:
                SeqIO.write(record, fh, "fasta")
            print(f"[✓] Reference saved → {REFERENCE_FASTA}  ({len(record.seq):,} bp)")
            return REFERENCE_FASTA
        except Exception as exc:
            print(f"[!] Attempt {attempt+1}/{retries} failed: {exc}")
            time.sleep(5)

    raise RuntimeError(f"Could not download reference genome after {retries} attempts.")


# ── Nextstrain open metadata ───────────────────────────────────────────────

def fetch_nextstrain_metadata(max_sequences: int = 10_000) -> Path:
    """
    Download Nextstrain open metadata TSV (zstd-compressed).
    Filters to *max_sequences* rows and saves as Parquet for fast loading.
    """
    parquet_path = OUTPUT_DIR / "nextstrain_metadata.parquet"

    if parquet_path.exists():
        print(f"[✓] Nextstrain metadata already present: {parquet_path}")
        return parquet_path

    try:
        import zstandard as zstd
    except ImportError:
        print("[!] zstandard not installed — falling back to sample data generator.")
        return _generate_sample_metadata(max_sequences)

    zst_path = RAW_DIR / "metadata.tsv.zst"
    if not zst_path.exists():
        print(f"[→] Downloading Nextstrain metadata (~500 MB, be patient) …")
        download_url(NEXTSTRAIN_METADATA_URL, zst_path, "nextstrain_metadata")

    print("[→] Decompressing and filtering metadata …")
    dctx = zstd.ZstdDecompressor()
    rows = []
    with open(zst_path, "rb") as fh:
        stream = dctx.stream_reader(fh)
        reader = io.TextIOWrapper(stream, encoding="utf-8")
        header = None
        n_cols = 0
        for i, line in enumerate(reader):
            if i == 0:
                header = line.strip().split("\t")
                n_cols = len(header)
                continue
            parts = line.strip().split("\t")
            # Pad or trim so every row has exactly n_cols fields
            if len(parts) < n_cols:
                parts += [""] * (n_cols - len(parts))
            elif len(parts) > n_cols:
                parts = parts[:n_cols]
            rows.append(parts)
            if len(rows) >= max_sequences:
                break

    df = pd.DataFrame(rows, columns=header)
    print(f"    Raw columns: {list(df.columns[:8])} ...")

    # Rename to standard internal column names
    col_map = {
        "Nextclade_pango":  "pangolin_lineage",
        "pango_lineage":    "pangolin_lineage",
        "Nextstrain_clade": "nextstrain_clade",
        "clade_nextstrain": "nextstrain_clade",
        "Nextclade_clade":  "nextstrain_clade",
        "clade_who":        "who_label",
    }
    for old_col, new_col in col_map.items():
        if old_col in df.columns and new_col not in df.columns:
            df = df.rename(columns={old_col: new_col})

    # Keep useful columns — accept whatever subset is present
    want = [
        "strain", "virus", "gisaid_epi_isl", "genbank_accession",
        "date", "region", "country", "division",
        "pangolin_lineage", "nextstrain_clade",
        "QC_overall_status",
        "substitutions", "deletions", "insertions",
        "aaSubstitutions", "aaDeletions",
    ]
    keep = [c for c in want if c in df.columns]
    df = df[keep]

    # Parse date — handle full dates (2023-01-15), month-only (2023-01), year-only (2023)
    def parse_flexible_date(d):
        if pd.isna(d) or not isinstance(d, str):
            return pd.NaT
        parts = d.strip().split("-")
        if len(parts) == 3:
            return pd.to_datetime(d, errors="coerce")
        elif len(parts) == 2:
            return pd.to_datetime(d + "-01", errors="coerce")
        elif len(parts) == 1 and len(d) == 4:
            return pd.to_datetime(d + "-01-01", errors="coerce")
        return pd.NaT
    df["date"] = df["date"].apply(parse_flexible_date)
    df = df.dropna(subset=["date"])
    df = df.sort_values("date")

    df.to_parquet(parquet_path, index=False)
    print(f"[✓] Metadata saved → {parquet_path}  ({len(df):,} rows)")
    return parquet_path


# ── Sample data generator (fallback when Nextstrain download not feasible) ─

def _generate_sample_metadata(n: int = 5000) -> Path:
    """
    Generate a realistic synthetic metadata table so the full pipeline
    runs even without internet / GISAID access.
    Uses real mutation patterns for each variant.
    """
    import random
    import numpy as np

    parquet_path = OUTPUT_DIR / "nextstrain_metadata.parquet"

    random.seed(42)
    np.random.seed(42)

    # Variant frequency schedule: (lineage, start_month, peak_month, end_month, peak_freq)
    variant_schedule = [
        ("Ancestral",      "2019-12", "2020-03", "2020-08", 0.90),
        ("B.1.1.7",        "2020-09", "2021-03", "2021-07", 0.85),
        ("B.1.351",        "2020-10", "2021-02", "2021-06", 0.30),
        ("P.1",            "2020-11", "2021-03", "2021-07", 0.25),
        ("B.1.617.2",      "2021-04", "2021-08", "2022-01", 0.95),
        ("B.1.1.529",      "2021-11", "2022-02", "2022-06", 0.90),
        ("BA.2",           "2021-12", "2022-04", "2022-09", 0.80),
        ("BA.4",           "2022-03", "2022-07", "2022-12", 0.40),
        ("BA.5",           "2022-03", "2022-08", "2023-01", 0.70),
        ("XBB.1.5",        "2022-10", "2023-02", "2023-08", 0.80),
        ("EG.5.1",         "2023-05", "2023-09", "2023-12", 0.60),
        ("JN.1",           "2023-08", "2024-01", "2024-06", 0.85),
    ]

    # Spike mutations per lineage
    spike_muts_by_lineage = {
        "Ancestral":   ["D614"],
        "B.1.1.7":     ["D614G", "N501Y", "P681H", "HV69-70del", "Y144del"],
        "B.1.351":     ["D614G", "K417N", "E484K", "N501Y"],
        "P.1":         ["D614G", "K417T", "E484K", "N501Y"],
        "B.1.617.2":   ["D614G", "L452R", "T478K", "P681R"],
        "B.1.1.529":   ["D614G", "K417N", "E484A", "N501Y", "P681H", "G339D",
                        "S371L", "Q493R", "G496S", "Q498R", "Y505H", "H655Y",
                        "N679K", "N764K", "D796Y"],
        "BA.2":        ["D614G", "K417N", "N501Y", "T478K", "G339D", "S371F",
                        "S373P", "Q493R", "Q498R", "Y505H"],
        "BA.4":        ["D614G", "K417N", "N501Y", "L452R", "F486V", "Y505H"],
        "BA.5":        ["D614G", "K417N", "N501Y", "L452R", "F486V", "Y505H"],
        "XBB.1.5":     ["D614G", "G339H", "R346T", "L368I", "V445P",
                        "N460K", "F486P", "F490S", "R493Q"],
        "EG.5.1":      ["D614G", "R346T", "N460K", "F456L", "F486P"],
        "JN.1":        ["D614G", "L455S", "F456L", "K417N", "N501Y",
                        "P681H", "N460K", "R346T"],
    }

    countries_by_region = {
        "Europe":   ["United Kingdom", "Germany", "France", "Denmark",
                     "Netherlands", "Sweden", "Italy", "Spain"],
        "Asia":     ["India", "Japan", "South Korea", "China",
                     "Singapore", "Thailand", "Israel"],
        "Americas": ["United States", "Brazil", "Canada", "Argentina",
                     "Mexico", "Colombia"],
        "Africa":   ["South Africa", "Kenya", "Nigeria", "Egypt"],
        "Oceania":  ["Australia", "New Zealand"],
    }

    all_dates = pd.date_range("2020-01-01", "2024-06-30", freq="D")

    records = []
    for i in range(n):
        date = random.choice(all_dates)
        ym   = date.strftime("%Y-%m")

        # Pick lineage by approximate schedule weights
        weights = []
        for lin, start, peak, end, max_f in variant_schedule:
            s = pd.Timestamp(start + "-01")
            p = pd.Timestamp(peak  + "-01")
            e = pd.Timestamp(end   + "-01")
            if date < s or date > e:
                weights.append(0.0)
            elif date <= p:
                frac = (date - s).days / max(1, (p - s).days)
                weights.append(max_f * frac)
            else:
                frac = 1.0 - (date - p).days / max(1, (e - p).days)
                weights.append(max_f * frac)

        w_arr = np.array(weights)
        total = w_arr.sum()
        if total == 0:
            w_arr = np.ones(len(weights))
        probs = w_arr / w_arr.sum()

        lineage = random.choices(
            [v[0] for v in variant_schedule], weights=probs
        )[0]

        region  = random.choices(
            list(countries_by_region.keys()),
            weights=[40, 25, 20, 10, 5]
        )[0]
        country = random.choice(countries_by_region[region])

        muts = spike_muts_by_lineage.get(lineage, [])
        # Add a few random background mutations
        bg = random.randint(0, 3)
        bg_muts = [f"{random.choice('ACDEFGHIKLMNPQRSTVWY')}{random.randint(1,1273)}{random.choice('ACDEFGHIKLMNPQRSTVWY')}"
                   for _ in range(bg)]

        aa_subs = ",".join(["S:" + m for m in muts] + ["S:" + m for m in bg_muts])

        records.append({
            "strain":           f"EPI_ISL_{1000000+i}",
            "date":             date,
            "region":           region,
            "country":          country,
            "pangolin_lineage": lineage,
            "nextstrain_clade": lineage,
            "aaSubstitutions":  aa_subs,
            "QC_overall_status": random.choices(["good","mediocre","bad"], weights=[85,10,5])[0],
        })

    df = pd.DataFrame(records)
    df = df[df["QC_overall_status"].isin(["good", "mediocre"])]
    df.to_parquet(parquet_path, index=False)
    print(f"[✓] Synthetic sample metadata saved → {parquet_path}  ({len(df):,} rows)")
    return parquet_path


# ── GISAID local FASTA (optional) ──────────────────────────────────────────

def process_gisaid_fasta(fasta_path: str, max_sequences: int = 2000) -> Path:
    """
    Parse a locally downloaded GISAID FASTA.
    Extracts metadata from sequence headers and saves a Parquet table.
    GISAID header format:
        >hCoV-19/Country/SampleID/Year|EPI_ISL_XXXXXX|2021-05-01
    """
    parquet_path = OUTPUT_DIR / "gisaid_metadata.parquet"
    records = []
    print(f"[→] Parsing GISAID FASTA: {fasta_path}")
    with open(fasta_path) as fh:
        for record in SeqIO.parse(fh, "fasta"):
            parts = record.description.split("|")
            if len(parts) < 3:
                continue
            name_parts = parts[0].split("/")
            country = name_parts[1] if len(name_parts) > 1 else "Unknown"
            epi     = parts[1].strip()
            date_s  = parts[2].strip()
            records.append({
                "strain":  record.id,
                "gisaid_epi_isl": epi,
                "country": country,
                "date":    pd.to_datetime(date_s, errors="coerce"),
                "sequence": str(record.seq),
            })
            if len(records) >= max_sequences:
                break

    df = pd.DataFrame(records).dropna(subset=["date"])
    df.to_parquet(parquet_path, index=False)
    print(f"[✓] GISAID metadata saved → {parquet_path}  ({len(df):,} rows)")
    return parquet_path


# ── CLI ────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Fetch COVID-19 genomic data")
    parser.add_argument("--source", choices=["nextstrain", "gisaid", "sample"],
                        default="sample",
                        help="Data source (default: sample — generates synthetic data)")
    parser.add_argument("--max-sequences", type=int, default=5000,
                        help="Maximum sequences to load (default: 5000)")
    parser.add_argument("--gisaid-fasta", type=str, default=None,
                        help="Path to locally downloaded GISAID sequences.fasta")
    args = parser.parse_args()

    fetch_reference_genome()

    if args.source == "gisaid":
        if not args.gisaid_fasta:
            parser.error("--gisaid-fasta required when --source gisaid")
        process_gisaid_fasta(args.gisaid_fasta, args.max_sequences)

    elif args.source == "nextstrain":
        fetch_nextstrain_metadata(args.max_sequences)

    else:  # sample
        print("[→] Generating synthetic sample data …")
        _generate_sample_metadata(args.max_sequences)

    print("\n[✓] Data fetch complete.")


if __name__ == "__main__":
    main()
