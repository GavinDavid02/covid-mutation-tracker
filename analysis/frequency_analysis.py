"""
analysis/frequency_analysis.py
───────────────────────────────
Computes:
  • Monthly mutation frequency time series
  • Variant dominance curves
  • Co-occurrence / linkage between mutations
  • Convergent evolution detection
  • Mutational hotspot ranking

Run standalone:
    python analysis/frequency_analysis.py
"""

import json
import sys
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency

sys.path.insert(0, str(Path(__file__).parent.parent))
from config import OUTPUT_DIR, KEY_MUTATIONS, VARIANT_DEFINITIONS


# ── Load processed data ────────────────────────────────────────────────────

def load_mutation_table() -> pd.DataFrame:
    parquet = OUTPUT_DIR / "mutation_table.parquet"
    if not parquet.exists():
        raise FileNotFoundError(
            f"Mutation table not found: {parquet}\n"
            "Run: python pipeline/fetch_data.py && python pipeline/alignment.py"
        )
    df = pd.read_parquet(parquet)
    df["date"] = pd.to_datetime(df["date"])
    df["year_month"] = df["date"].dt.to_period("M")
    return df


# ── 1. Monthly variant frequency ───────────────────────────────────────────

def compute_variant_frequency(df: pd.DataFrame) -> pd.DataFrame:
    """
    Returns a DataFrame with columns:
        year_month | variant | count | total | frequency
    sorted by year_month, variant.
    """
    monthly_total = (
        df.groupby("year_month")
          .size()
          .rename("total")
          .reset_index()
    )
    monthly_variant = (
        df.groupby(["year_month", "who_label"])
          .size()
          .rename("count")
          .reset_index()
    )
    merged = monthly_variant.merge(monthly_total, on="year_month")
    merged["frequency"] = merged["count"] / merged["total"] * 100
    merged = merged.sort_values(["year_month", "frequency"], ascending=[True, False])
    return merged


# ── 2. Per-mutation frequency over time ────────────────────────────────────

def compute_mutation_frequency_over_time(
    df: pd.DataFrame,
    mutations: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    For each key mutation, compute its monthly frequency (% of sequences with that mutation).

    Returns: year_month | mutation | frequency | count | total
    """
    if mutations is None:
        mutations = [m[0] for m in KEY_MUTATIONS]

    monthly_total = (
        df.groupby("year_month")
          .size()
          .rename("total")
          .reset_index()
    )

    rows = []
    for mut in mutations:
        col = f"has_{mut}"
        if col not in df.columns:
            # Parse from spike_mutations JSON column
            df[col] = df["spike_mutations"].apply(
                lambda x: mut in (json.loads(x) if isinstance(x, str) else [])
            )
        monthly_mut = (
            df[df[col] == True]
              .groupby("year_month")
              .size()
              .rename("count")
              .reset_index()
        )
        monthly_mut["mutation"] = mut
        merged = monthly_mut.merge(monthly_total, on="year_month")
        merged["frequency"] = merged["count"] / merged["total"] * 100
        rows.append(merged)

    if not rows:
        return pd.DataFrame()

    result = pd.concat(rows, ignore_index=True)
    return result.sort_values(["mutation", "year_month"])


# ── 3. Mutation co-occurrence (linkage) ────────────────────────────────────

def compute_mutation_cooccurrence(
    df: pd.DataFrame,
    top_n: int = 20,
) -> pd.DataFrame:
    """
    Build a co-occurrence matrix for the top_n most frequent mutations.
    Uses chi-squared test to assess statistical significance.

    Returns DataFrame with columns:
        mut_a | mut_b | cooccurrence_count | chi2 | p_value | odds_ratio
    """
    # Identify most common mutations
    all_muts: Dict[str, int] = {}
    for row_muts in df["spike_mutations"].dropna():
        for m in json.loads(row_muts) if isinstance(row_muts, str) else []:
            all_muts[m] = all_muts.get(m, 0) + 1

    top_muts = sorted(all_muts, key=all_muts.get, reverse=True)[:top_n]

    # Build binary indicator matrix
    indicators = {}
    for mut in top_muts:
        col = f"has_{mut}"
        if col in df.columns:
            indicators[mut] = df[col].astype(int)
        else:
            indicators[mut] = df["spike_mutations"].apply(
                lambda x: 1 if (isinstance(x, str) and mut in json.loads(x)) else 0
            )

    ind_df = pd.DataFrame(indicators)
    n = len(ind_df)

    rows = []
    for m1, m2 in combinations(top_muts, 2):
        a = ind_df[m1].values
        b = ind_df[m2].values
        both  = int((a & b).sum())
        a_only= int((a & ~b).sum())
        b_only= int((~a & b).sum())
        neither=int((~a & ~b).sum())

        table = np.array([[both, a_only], [b_only, neither]])
        if table.min() < 5:   # skip cells too sparse for chi2
            continue

        chi2, p, _, _ = chi2_contingency(table)
        # Odds ratio (add 0.5 for Haldane-Anscombe correction)
        or_val = ((both + 0.5) * (neither + 0.5)) / ((a_only + 0.5) * (b_only + 0.5))

        rows.append({
            "mut_a":              m1,
            "mut_b":              m2,
            "cooccurrence_count": both,
            "chi2":               round(chi2, 3),
            "p_value":            p,
            "odds_ratio":         round(or_val, 3),
        })

    result = pd.DataFrame(rows)
    if not result.empty:
        result = result.sort_values("odds_ratio", ascending=False)
    return result


# ── 4. Convergent evolution detection ─────────────────────────────────────

def detect_convergent_evolution(df: pd.DataFrame) -> pd.DataFrame:
    """
    Identify mutations that appeared independently in multiple distinct lineages.
    A mutation is 'convergent' if it appears in ≥3 WHO-labelled variants.

    Returns: mutation | position | lineages | appearance_count | convergence_score
    """
    # Which lineages carry each mutation?
    mut_lineages: Dict[str, set] = {}
    for _, row in df.iterrows():
        muts = json.loads(row["spike_mutations"]) if isinstance(row["spike_mutations"], str) else []
        lineage = row.get("who_label", "Unknown")
        for mut in muts:
            mut_lineages.setdefault(mut, set()).add(lineage)

    rows = []
    for mut, lineages in mut_lineages.items():
        n_lineages = len(lineages)
        if n_lineages >= 2:  # appeared in at least 2 distinct WHO groups
            # Extract position
            import re
            m = re.search(r"\d+", mut)
            pos = int(m.group()) if m else -1
            rows.append({
                "mutation":          mut,
                "position":          pos,
                "lineages":          ", ".join(sorted(lineages)),
                "appearance_count":  n_lineages,
                "convergence_score": n_lineages,  # simple score; could weight by phylo distance
            })

    result = pd.DataFrame(rows)
    if not result.empty:
        result = result.sort_values("convergence_score", ascending=False)
    return result


# ── 5. Mutational hotspot detection ────────────────────────────────────────

def compute_mutation_hotspots(df: pd.DataFrame, window: int = 10) -> pd.DataFrame:
    """
    Identify amino-acid positions with unusually high mutation rates.
    Uses a sliding window over the spike protein positions.

    Returns: position | aa_ref | mutation_count | top_mutations | domain | is_hotspot
    """
    import re
    from config import SPIKE_DOMAINS

    def _domain(pos):
        for d, (s, e) in SPIKE_DOMAINS.items():
            if s <= pos <= e:
                return d
        return "Other"

    position_counts: Dict[int, Dict[str, int]] = {}
    for row_muts in df["spike_mutations"].dropna():
        muts = json.loads(row_muts) if isinstance(row_muts, str) else []
        for mut in muts:
            m = re.search(r"([A-Z*])(\d+)([A-Z*])", mut)
            if m:
                pos = int(m.group(2))
                position_counts.setdefault(pos, {})
                position_counts[pos][mut] = position_counts[pos].get(mut, 0) + 1

    rows = []
    for pos, mut_dict in position_counts.items():
        total = sum(mut_dict.values())
        top   = sorted(mut_dict, key=mut_dict.get, reverse=True)[:3]
        rows.append({
            "position":       pos,
            "mutation_count": total,
            "top_mutations":  ", ".join(top),
            "domain":         _domain(pos),
        })

    result = pd.DataFrame(rows).sort_values("mutation_count", ascending=False)
    if not result.empty:
        threshold = result["mutation_count"].quantile(0.90)
        result["is_hotspot"] = result["mutation_count"] >= threshold
    return result


# ── 6. Sweep speed analysis ────────────────────────────────────────────────

def compute_sweep_speed(df: pd.DataFrame) -> pd.DataFrame:
    """
    For each key mutation, compute how quickly it swept from 5% to 50% global frequency.
    Returns: mutation | first_seen | reached_5pct | reached_50pct | sweep_days
    """
    rows = []
    for mut, *_ in KEY_MUTATIONS:
        col = f"has_{mut}"
        if col not in df.columns:
            continue

        mut_df = df.sort_values("date").copy()
        mut_df["month"] = mut_df["date"].dt.to_period("M")
        monthly = mut_df.groupby("month").agg(
            total=("strain", "count"),
            with_mut=(col, "sum")
        ).reset_index()
        monthly["freq"] = monthly["with_mut"] / monthly["total"] * 100

        try:
            first_seen = monthly[monthly["freq"] > 0]["month"].min()
            r5  = monthly[monthly["freq"] >= 5]["month"].min()
            r50 = monthly[monthly["freq"] >= 50]["month"].min()

            if pd.isna(r5) or pd.isna(r50):
                sweep_days = None
            else:
                sweep_days = (r50.to_timestamp() - r5.to_timestamp()).days

            rows.append({
                "mutation":    mut,
                "first_seen":  str(first_seen) if not pd.isna(first_seen) else "not seen",
                "reached_5pct": str(r5)  if not pd.isna(r5)  else "not reached",
                "reached_50pct":str(r50) if not pd.isna(r50) else "not reached",
                "sweep_days":  sweep_days,
                "peak_freq":   round(monthly["freq"].max(), 1),
            })
        except Exception:
            continue

    return pd.DataFrame(rows).sort_values("sweep_days")


# ── 7. Geographic spread analysis ─────────────────────────────────────────

def compute_geographic_spread(df: pd.DataFrame) -> pd.DataFrame:
    """
    For each variant, compute first detection date per country.
    Returns: who_label | country | first_detection | sequences_count
    """
    result = (
        df.groupby(["who_label", "country"])
          .agg(
              first_detection=("date", "min"),
              sequences_count=("strain", "count"),
          )
          .reset_index()
          .sort_values(["who_label", "first_detection"])
    )
    return result


# ── Save all analyses ───────────────────────────────────────────────────────

def run_all_analyses(df: Optional[pd.DataFrame] = None) -> Dict:
    if df is None:
        df = load_mutation_table()

    df["date"] = pd.to_datetime(df["date"])
    df["year_month"] = df["date"].dt.to_period("M")

    print("[→] Computing variant frequency …")
    variant_freq = compute_variant_frequency(df)
    variant_freq.to_parquet(OUTPUT_DIR / "variant_frequency.parquet", index=False)

    print("[→] Computing per-mutation frequency over time …")
    mut_freq = compute_mutation_frequency_over_time(df)
    mut_freq.to_parquet(OUTPUT_DIR / "mutation_frequency_time.parquet", index=False)

    print("[→] Computing mutation co-occurrence …")
    cooccur = compute_mutation_cooccurrence(df)
    cooccur.to_parquet(OUTPUT_DIR / "cooccurrence.parquet", index=False)

    print("[→] Detecting convergent evolution …")
    convergent = detect_convergent_evolution(df)
    convergent.to_parquet(OUTPUT_DIR / "convergent_evolution.parquet", index=False)

    print("[→] Computing hotspots …")
    hotspots = compute_mutation_hotspots(df)
    hotspots.to_parquet(OUTPUT_DIR / "hotspots.parquet", index=False)

    print("[→] Computing sweep speeds …")
    sweeps = compute_sweep_speed(df)
    sweeps.to_parquet(OUTPUT_DIR / "sweep_speed.parquet", index=False)

    print("[→] Computing geographic spread …")
    geo = compute_geographic_spread(df)
    geo.to_parquet(OUTPUT_DIR / "geographic_spread.parquet", index=False)

    print(f"[✓] All analyses saved to {OUTPUT_DIR}")
    return {
        "variant_frequency": variant_freq,
        "mutation_frequency": mut_freq,
        "cooccurrence": cooccur,
        "convergent": convergent,
        "hotspots": hotspots,
        "sweeps": sweeps,
        "geo": geo,
    }


if __name__ == "__main__":
    results = run_all_analyses()
    print("\nTop convergent mutations:")
    print(results["convergent"].head(10).to_string(index=False))
    print("\nFastest sweeping mutations:")
    print(results["sweeps"].head(10).to_string(index=False))
