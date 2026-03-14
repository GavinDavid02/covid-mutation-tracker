"""
analysis/generate_plots.py
──────────────────────────
Standalone script that runs the full analysis and exports every figure
as a self-contained HTML file (interactive Plotly) and a static PNG.

Run:
    python analysis/generate_plots.py
    python analysis/generate_plots.py --format html     # HTML only (no kaleido needed)
    python analysis/generate_plots.py --format png      # PNG only (requires kaleido)
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))


def export_figure(fig, name: str, output_dir: Path, fmt: str):
    output_dir.mkdir(parents=True, exist_ok=True)
    if fmt in ("html", "both"):
        path = output_dir / f"{name}.html"
        fig.write_html(str(path), include_plotlyjs="cdn")
        print(f"  [✓] {path.name}")
    if fmt in ("png", "both"):
        try:
            path = output_dir / f"{name}.png"
            fig.write_image(str(path), width=1200, height=600, scale=2)
            print(f"  [✓] {path.name}")
        except Exception as e:
            print(f"  [!] PNG export failed (install kaleido): {e}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--format", choices=["html", "png", "both"], default="html")
    args = parser.parse_args()

    from config import OUTPUT_DIR
    from analysis.frequency_analysis import (
        load_mutation_table,
        compute_variant_frequency,
        compute_mutation_frequency_over_time,
        detect_convergent_evolution,
        compute_mutation_hotspots,
        compute_sweep_speed,
        compute_geographic_spread,
        compute_mutation_cooccurrence,
    )
    from visualization.plots import (
        timeline_plot,
        geographic_plotly,
        spike_protein_map,
        spike_rbd_detail,
        phylogenetic_tree_plotly,
        cooccurrence_heatmap,
        sweep_speed_chart,
    )

    out_dir = OUTPUT_DIR / "figures"
    print(f"[→] Loading mutation table …")
    df = load_mutation_table()
    print(f"    {len(df):,} sequences loaded")

    print(f"\n[→] Running analyses …")
    vf      = compute_variant_frequency(df)
    mf      = compute_mutation_frequency_over_time(df)
    conv    = detect_convergent_evolution(df)
    hs      = compute_mutation_hotspots(df)
    sweeps  = compute_sweep_speed(df)
    geo     = compute_geographic_spread(df)
    cooccur = compute_mutation_cooccurrence(df)

    print(f"\n[→] Exporting figures → {out_dir}")

    figs = [
        (timeline_plot(vf, mode="variant"),                          "01_variant_dominance"),
        (timeline_plot(vf, mf, mode="mutation"),                     "02_mutation_frequency"),
        (geographic_plotly(geo),                                     "03_geographic_spread"),
        (spike_protein_map(hs),                                      "04_spike_protein_map"),
        (spike_rbd_detail(hs),                                       "05_rbd_zoom"),
        (phylogenetic_tree_plotly(),                                  "06_phylogenetic_tree"),
        (cooccurrence_heatmap(cooccur),                              "07_cooccurrence_heatmap"),
        (sweep_speed_chart(sweeps),                                  "08_sweep_speed"),
    ]

    for fig, name in figs:
        export_figure(fig, name, out_dir, args.format)

    print(f"\n[✓] {len(figs)} figures exported to {out_dir}")


if __name__ == "__main__":
    main()
