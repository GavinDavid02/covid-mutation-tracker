"""
pipeline/alignment.py
─────────────────────
Pairwise-align SARS-CoV-2 spike sequences to the Wuhan-Hu-1 reference and
call amino-acid mutations.

Two alignment back-ends are supported:
  • Biopython PairwiseAligner  (pure Python, no external tools, always available)
  • minimap2                   (much faster; used automatically if installed)

Outputs a per-sample mutation table:
    strain | date | country | lineage | spike_mutations | rbd_mutations | rbd_count | ...

Run:
    python pipeline/alignment.py
    python pipeline/alignment.py --backend minimap2 --workers 4
"""

import argparse
import json
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import List, Optional, Tuple

import pandas as pd
import numpy as np
from Bio import SeqIO, pairwise2
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).parent.parent))
from config import (
    REFERENCE_FASTA, OUTPUT_DIR, RAW_DIR,
    SPIKE_START, SPIKE_END, SPIKE_DOMAINS,
    KEY_MUTATIONS,
)

# ── Reference helpers ──────────────────────────────────────────────────────

def load_reference() -> str:
    """Return the full reference genome sequence as a string."""
    if not REFERENCE_FASTA.exists():
        raise FileNotFoundError(
            f"Reference not found: {REFERENCE_FASTA}\n"
            "Run: python pipeline/fetch_data.py"
        )
    record = next(SeqIO.parse(REFERENCE_FASTA, "fasta"))
    return str(record.seq).upper()


def extract_spike_aa_reference(ref_nt: str) -> str:
    """
    Translate the spike ORF from the Wuhan-Hu-1 reference.
    Uses 0-based Python slicing; coordinates in config are 1-based.
    """
    spike_nt = ref_nt[SPIKE_START - 1: SPIKE_END]
    codon_count = len(spike_nt) // 3
    aa_seq = ""
    for i in range(codon_count):
        codon = spike_nt[i*3: i*3+3]
        codon_seq = Seq(codon)
        aa_seq += str(codon_seq.translate())
    return aa_seq


# ── Mutation calling from Nextclade / Nextstrain aaSubstitutions ───────────

def parse_aa_substitutions(aa_subs_str: str, gene: str = "S") -> List[str]:
    """
    Parse aaSubstitutions field from Nextstrain metadata.
    Format: "S:D614G,S:N501Y,ORF1a:T4715I"
    Returns list of spike substitutions like ["D614G", "N501Y"].
    """
    if not aa_subs_str or pd.isna(aa_subs_str):
        return []
    mutations = []
    for part in str(aa_subs_str).split(","):
        part = part.strip()
        if part.startswith(gene + ":"):
            mut = part[len(gene)+1:]
            mutations.append(mut)
    return mutations


def parse_aa_deletions(aa_del_str: str, gene: str = "S") -> List[str]:
    """
    Parse aaDeletions field.
    Format: "S:HV69-70,ORF1a:3675-3677"
    Returns list of deletion strings like ["HV69-70del"].
    """
    if not aa_del_str or pd.isna(aa_del_str):
        return []
    deletions = []
    for part in str(aa_del_str).split(","):
        part = part.strip()
        if part.startswith(gene + ":"):
            pos = part[len(gene)+1:]
            deletions.append(pos + "del")
    return deletions


# ── Biopython pairwise aligner (fallback) ─────────────────────────────────

def align_spike_biopython(query_nt: str, ref_nt: str) -> List[str]:
    """
    Align a query spike nucleotide sequence against the reference,
    translate both, and return list of amino-acid substitutions.
    """
    ref_spike = ref_nt[SPIKE_START - 1: SPIKE_END]
    # Trim query to spike region (heuristic: first 3kb of query assumed to be spike)
    query_spike = query_nt[:len(ref_spike)]

    aligner = PairwiseAligner()
    aligner.mode             = "global"
    aligner.match_score      = 2
    aligner.mismatch_score   = -1
    aligner.open_gap_score   = -5
    aligner.extend_gap_score = -0.5

    try:
        alignments = aligner.align(ref_spike, query_spike)
        alignment  = next(iter(alignments))
        ref_aligned   = alignment[0]
        query_aligned = alignment[1]
    except Exception:
        return []

    # Walk aligned sequences, translate codons, find substitutions
    ref_pos   = 0   # nucleotide position in ref
    query_pos = 0
    mutations = []

    ref_codons   = []
    query_codons = []
    for r_nt, q_nt in zip(ref_aligned, query_aligned):
        if r_nt != "-":
            ref_codons.append(r_nt)
        if q_nt != "-":
            query_codons.append(q_nt)

    ref_aa   = Seq("".join(ref_codons)).translate()
    query_aa = Seq("".join(query_codons)).translate()

    for i, (r, q) in enumerate(zip(ref_aa, query_aa)):
        if r != q and q != "X":
            mutations.append(f"{r}{i+1}{q}")

    return mutations


# ── minimap2 wrapper ────────────────────────────────────────────────────────

def _minimap2_available() -> bool:
    try:
        result = subprocess.run(["minimap2", "--version"],
                                capture_output=True, text=True)
        return result.returncode == 0
    except FileNotFoundError:
        return False


def align_spike_minimap2(sequences: List[Tuple[str, str]], ref_fasta: Path) -> dict:
    """
    Align a batch of (strain, sequence) tuples using minimap2.
    Returns dict: strain -> list of mutations.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        query_fasta = tmpdir / "query.fasta"
        with open(query_fasta, "w") as fh:
            for strain, seq in sequences:
                fh.write(f">{strain}\n{seq}\n")

        sam_out = tmpdir / "aln.sam"
        cmd = [
            "minimap2", "-a", "--cs",
            "-x", "map-ont",   # works for ~30 kb coronavirus genomes
            str(ref_fasta), str(query_fasta)
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"minimap2 failed: {result.stderr[:500]}")

        with open(sam_out, "w") as fh:
            fh.write(result.stdout)

        # Parse CS strings from SAM to extract mutations (simplified)
        mutations_by_strain = {}
        for line in result.stdout.splitlines():
            if line.startswith("@"):
                continue
            fields = line.split("\t")
            if len(fields) < 11:
                continue
            qname = fields[0]
            # Extract cs tag
            cs_tag = ""
            for f in fields[11:]:
                if f.startswith("cs:Z:"):
                    cs_tag = f[5:]
                    break
            mutations_by_strain[qname] = _parse_cs_to_mutations(cs_tag)

        return mutations_by_strain


def _parse_cs_to_mutations(cs: str) -> List[str]:
    """
    Parse minimap2 cs string to extract substitution positions.
    cs format: :42*at+tcc:10 → match 42, sub A→T, ins TCC, match 10
    """
    import re
    mutations = []
    pos = 0
    for token in re.findall(r"[=:*+-][A-Za-z0-9]+", cs):
        op = token[0]
        val = token[1:]
        if op in (":", "="):
            pos += len(val)
        elif op == "*":
            ref_aa = val[0].upper()
            alt_aa = val[1].upper()
            mutations.append(f"{ref_aa}{pos+1}{alt_aa}")
            pos += 1
        elif op == "+":
            pass  # insertion
        elif op == "-":
            pos += len(val)
    return mutations


# ── Main processing: metadata → mutation table ─────────────────────────────

def build_mutation_table(
    metadata_parquet: Path,
    backend: str = "metadata",
    workers: int = 2,
) -> pd.DataFrame:
    """
    Build the core per-sample mutation table from metadata.

    backend = "metadata"  → parse aaSubstitutions column (fast, no alignment)
    backend = "biopython" → pairwise-align spike sequences (slow, needs FASTA)
    backend = "minimap2"  → minimap2 alignment (fast, needs binary + FASTA)
    """
    print(f"[→] Building mutation table from {metadata_parquet} …")
    df = pd.read_parquet(metadata_parquet)

    records = []
    for _, row in tqdm(df.iterrows(), total=len(df), desc="Parsing mutations"):
        lineage = str(row.get("pangolin_lineage", row.get("nextstrain_clade", "Unknown")))
        aa_subs = row.get("aaSubstitutions", "")
        aa_dels = row.get("aaDeletions", "")

        spike_muts = parse_aa_substitutions(aa_subs, gene="S")
        spike_dels = parse_aa_deletions(aa_dels, gene="S")
        all_spike  = spike_muts + spike_dels

        # Annotate each mutation with its domain
        mut_with_domains = []
        for mut in all_spike:
            pos = _extract_position(mut)
            domain = _position_to_domain(pos) if pos else "Unknown"
            mut_with_domains.append({"mutation": mut, "position": pos, "domain": domain})

        # RBD mutations specifically
        rbd_muts = [m["mutation"] for m in mut_with_domains
                    if m["domain"] in ("RBD", "RBM")]

        # Check for key known mutations
        key_mut_flags = {}
        for km, pos, cat, _ in KEY_MUTATIONS:
            key_mut_flags[f"has_{km}"] = km in all_spike

        records.append({
            "strain":          row.get("strain", f"seq_{_}"),
            "date":            row.get("date"),
            "region":          row.get("region", "Unknown"),
            "country":         row.get("country", "Unknown"),
            "lineage":         lineage,
            "who_label":       _lineage_to_who(lineage),
            "spike_mutations": json.dumps(all_spike),
            "spike_mut_count": len(all_spike),
            "rbd_mutations":   json.dumps(rbd_muts),
            "rbd_mut_count":   len(rbd_muts),
            **key_mut_flags,
        })

    result_df = pd.DataFrame(records)
    result_df["date"] = pd.to_datetime(result_df["date"])
    result_df = result_df.sort_values("date").reset_index(drop=True)

    out_path = OUTPUT_DIR / "mutation_table.parquet"
    result_df.to_parquet(out_path, index=False)
    print(f"[✓] Mutation table saved → {out_path}  ({len(result_df):,} rows)")
    return result_df


# ── Helpers ────────────────────────────────────────────────────────────────

def _extract_position(mutation: str) -> Optional[int]:
    """Extract numeric position from a mutation string like 'D614G' → 614."""
    import re
    m = re.search(r"\d+", mutation)
    return int(m.group()) if m else None


def _position_to_domain(pos: int) -> str:
    """Return the spike domain for a given amino-acid position."""
    for domain, (start, end) in SPIKE_DOMAINS.items():
        if start <= pos <= end:
            return domain
    return "Other"


def _lineage_to_who(lineage: str) -> str:
    """Map a Pango lineage string to its WHO label."""
    mapping = {
        "B.1.1.7":  "Alpha",
        "B.1.351":  "Beta",
        "P.1":      "Gamma",
        "B.1.617.2":"Delta",
        "B.1.1.529":"Omicron",
        "BA.":      "Omicron",
        "BQ.":      "Omicron",
        "XBB":      "Omicron",
        "EG.":      "Omicron",
        "JN.":      "Omicron",
        "KP.":      "Omicron",
    }
    lin = str(lineage).strip()
    for k, v in mapping.items():
        if lin.startswith(k):
            return v
    if lin in ("Ancestral", "Unknown", "nan", ""):
        return "Ancestral"
    return "Other"


# ── CLI ────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Build mutation table from metadata")
    parser.add_argument("--backend",
                        choices=["metadata", "biopython", "minimap2"],
                        default="metadata")
    parser.add_argument("--workers", type=int, default=2)
    args = parser.parse_args()

    parquet = OUTPUT_DIR / "nextstrain_metadata.parquet"
    if not parquet.exists():
        print("[!] Metadata not found — run pipeline/fetch_data.py first.")
        sys.exit(1)

    build_mutation_table(parquet, backend=args.backend, workers=args.workers)


if __name__ == "__main__":
    main()
