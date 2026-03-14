"""
pipeline/run_pipeline.py
─────────────────────────
End-to-end pipeline runner.

Usage:
    # Quick start with synthetic data (no internet required):
    python pipeline/run_pipeline.py --mode sample

    # Real data from Nextstrain (requires internet):
    python pipeline/run_pipeline.py --mode nextstrain --max-sequences 50000

    # Real data from GISAID (requires local FASTA file):
    python pipeline/run_pipeline.py --mode gisaid --gisaid-fasta /path/to/sequences.fasta

    # Skip steps you've already run:
    python pipeline/run_pipeline.py --mode nextstrain --skip-fetch --skip-align
"""

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))


def print_banner():
    print("""
╔══════════════════════════════════════════════════════╗
║   SARS-CoV-2 Mutation Tracker — Pipeline Runner     ║
║   Wuhan reference · Nextstrain · GISAID              ║
╚══════════════════════════════════════════════════════╝
""")


def step(n, total, name):
    print(f"\n[Step {n}/{total}] {name}")
    print("─" * 50)


def main():
    print_banner()

    parser = argparse.ArgumentParser(
        description="Run the full COVID-19 mutation tracking pipeline"
    )
    parser.add_argument(
        "--mode",
        choices=["sample", "nextstrain", "gisaid"],
        default="sample",
        help="Data source (default: sample)"
    )
    parser.add_argument(
        "--max-sequences", type=int, default=10_000,
        help="Max sequences to process (default: 10000)"
    )
    parser.add_argument(
        "--gisaid-fasta", type=str, default=None,
        help="Path to GISAID sequences.fasta (required for gisaid mode)"
    )
    parser.add_argument(
        "--alignment-backend",
        choices=["metadata", "biopython", "minimap2"],
        default="metadata",
        help="Alignment backend (default: metadata)"
    )
    parser.add_argument(
        "--workers", type=int, default=2,
        help="Worker processes for alignment (default: 2)"
    )
    parser.add_argument("--skip-fetch",  action="store_true")
    parser.add_argument("--skip-align",  action="store_true")
    parser.add_argument("--skip-analyze",action="store_true")

    args = parser.parse_args()

    total_steps = 3 - sum([args.skip_fetch, args.skip_align, args.skip_analyze])
    current = 0
    t0 = time.time()

    from config import OUTPUT_DIR

    # ── Step 1: Fetch data ─────────────────────────────────────────────────
    if not args.skip_fetch:
        current += 1
        step(current, total_steps, "Fetching data")
        from pipeline.fetch_data import (
            fetch_reference_genome,
            fetch_nextstrain_metadata,
            process_gisaid_fasta,
            _generate_sample_metadata,
        )
        fetch_reference_genome()

        if args.mode == "sample":
            meta_path = _generate_sample_metadata(args.max_sequences)
        elif args.mode == "nextstrain":
            meta_path = fetch_nextstrain_metadata(args.max_sequences)
        elif args.mode == "gisaid":
            if not args.gisaid_fasta:
                print("[!] --gisaid-fasta is required for gisaid mode")
                sys.exit(1)
            meta_path = process_gisaid_fasta(args.gisaid_fasta, args.max_sequences)
    else:
        print("[skip] Data fetch")
        meta_path = OUTPUT_DIR / "nextstrain_metadata.parquet"

    # ── Step 2: Alignment & mutation calling ──────────────────────────────
    if not args.skip_align:
        current += 1
        step(current, total_steps, "Alignment & mutation calling")
        from pipeline.alignment import build_mutation_table
        df = build_mutation_table(
            meta_path,
            backend=args.alignment_backend,
            workers=args.workers,
        )
    else:
        print("[skip] Alignment")
        import pandas as pd
        mut_path = OUTPUT_DIR / "mutation_table.parquet"
        if not mut_path.exists():
            print("[!] mutation_table.parquet not found — cannot skip alignment")
            sys.exit(1)
        df = pd.read_parquet(mut_path)

    # ── Step 3: Frequency & statistical analysis ──────────────────────────
    if not args.skip_analyze:
        current += 1
        step(current, total_steps, "Frequency & statistical analysis")
        from analysis.frequency_analysis import run_all_analyses
        results = run_all_analyses(df)

        print("\n── Summary ──────────────────────────────────────────")
        print(f"Sequences processed:    {len(df):,}")
        print(f"Unique WHO variants:    {df['who_label'].nunique()}")
        print(f"Countries represented:  {df['country'].nunique()}")

        conv = results.get("convergent", None)
        if conv is not None and not conv.empty:
            top3 = conv.nlargest(3, "appearance_count")["mutation"].tolist()
            print(f"Top convergent muts:    {', '.join(top3)}")

        sweeps = results.get("sweeps", None)
        if sweeps is not None and not sweeps.empty:
            fastest = sweeps.dropna(subset=["sweep_days"]).nsmallest(1, "sweep_days")
            if not fastest.empty:
                r = fastest.iloc[0]
                print(f"Fastest sweep:          {r['mutation']} ({int(r['sweep_days'])} days)")

    else:
        print("[skip] Analysis")

    elapsed = time.time() - t0
    print(f"\n✓ Pipeline complete in {elapsed:.1f}s")
    print(f"  Processed data:  {OUTPUT_DIR}")
    print(f"\n  Launch dashboard:")
    print(f"    streamlit run app.py")


if __name__ == "__main__":
    main()
