"""
visualization/plots.py
──────────────────────
All four visualisation modules for the COVID-19 Mutation Tracker:

  1. timeline_plot()      — mutation / variant frequency curves over time
  2. geographic_map()     — folium choropleth + first-detection markers
  3. spike_protein_map()  — plotly linear map of spike with mutation hotspots
  4. phylogenetic_tree()  — ETE3 or plotly-based variant tree

Each function returns a Plotly figure (or a folium Map for geographic_map).
"""

import json
import sys
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

sys.path.insert(0, str(Path(__file__).parent.parent))
from config import (
    SPIKE_DOMAINS, KEY_MUTATIONS, VARIANT_DEFINITIONS,
    FUNCTIONAL_COLOURS, OUTPUT_DIR,
)


# ── Colour helpers ─────────────────────────────────────────────────────────

WHO_COLOURS = {v["who"]: v["colour"] for v in VARIANT_DEFINITIONS.values()}
WHO_COLOURS.update({
    "Ancestral": "#888780",
    "Other":     "#B4B2A9",
})

DOMAIN_COLOURS = {
    "Signal peptide": "#D3D1C7",
    "NTD":            "#85B7EB",
    "RBD":            "#F09595",
    "RBM":            "#E24B4A",
    "Central domain": "#B5D4F4",
    "Furin cleavage": "#EF9F27",
    "Fusion peptide": "#97C459",
    "HR1":            "#FAC775",
    "HR2":            "#9FE1CB",
    "TM":             "#C0DD97",
    "Cytoplasmic":    "#D3D1C7",
    "Other":          "#888780",
}


# ── 1. Timeline plot ───────────────────────────────────────────────────────

def timeline_plot(
    variant_freq: pd.DataFrame,
    mutation_freq: Optional[pd.DataFrame] = None,
    mode: str = "variant",       # "variant" | "mutation"
    selected_mutations: Optional[List[str]] = None,
    region: Optional[str] = None,
) -> go.Figure:
    """
    mode="variant"  → stacked area chart of WHO-labelled variant dominance
    mode="mutation" → line chart of individual mutation frequencies
    """
    fig = go.Figure()

    if mode == "variant":
        vf = variant_freq.copy()
        vf["year_month_ts"] = vf["year_month"].dt.to_timestamp()

        # Pivot to wide
        pivot = vf.pivot_table(
            index="year_month_ts",
            columns="who_label",
            values="frequency",
            aggfunc="sum",
        ).fillna(0)

        # Order variants chronologically by their first appearance
        order = ["Ancestral", "Alpha", "Beta", "Gamma", "Delta", "Omicron", "Other"]
        cols = [c for c in order if c in pivot.columns] + \
               [c for c in pivot.columns if c not in order]

        for col in cols:
            fig.add_trace(go.Scatter(
                x=pivot.index,
                y=pivot[col].round(1),
                name=col,
                mode="lines",
                line=dict(width=2, color=WHO_COLOURS.get(col, "#888")),
                fill="tonexty" if col != cols[0] else "tozeroy",
                stackgroup="one",
                hovertemplate=f"<b>{col}</b><br>%{{x|%b %Y}}<br>Frequency: %{{y:.1f}}%<extra></extra>",
            ))

        fig.update_layout(
            title="Global SARS-CoV-2 variant dominance over time",
            xaxis_title="Date",
            yaxis_title="% of sequences",
            yaxis=dict(range=[0, 100]),
            hovermode="x unified",
            legend=dict(orientation="h", y=-0.2),
        )

    elif mode == "mutation":
        if mutation_freq is None or mutation_freq.empty:
            return _empty_fig("No mutation frequency data available")

        mf = mutation_freq.copy()
        mf["year_month_ts"] = mf["year_month"].dt.to_timestamp()

        if selected_mutations:
            mf = mf[mf["mutation"].isin(selected_mutations)]

        for mut, grp in mf.groupby("mutation"):
            cat  = next((m[2] for m in KEY_MUTATIONS if m[0] == mut), "Other")
            col  = FUNCTIONAL_COLOURS.get(cat, "#888780")
            fig.add_trace(go.Scatter(
                x=grp["year_month_ts"],
                y=grp["frequency"].round(1),
                name=mut,
                mode="lines+markers",
                line=dict(width=2, color=col),
                marker=dict(size=4),
                hovertemplate=f"<b>{mut}</b><br>%{{x|%b %Y}}<br>Frequency: %{{y:.1f}}%<extra></extra>",
            ))

        fig.update_layout(
            title="Spike mutation frequency over time",
            xaxis_title="Date",
            yaxis_title="% of sequences",
            yaxis=dict(range=[0, 105]),
            hovermode="x unified",
            legend=dict(orientation="h", y=-0.25, font=dict(size=11)),
        )

    _apply_theme(fig)
    return fig


# ── 2. Geographic map ──────────────────────────────────────────────────────

def geographic_map(geo_df: pd.DataFrame) -> "folium.Map":
    """
    Build a folium choropleth + circle-marker map showing:
      • Country fill = number of sequences
      • Circle marker = first-detection event per variant
    """
    try:
        import folium
        from folium.plugins import MarkerCluster
    except ImportError:
        return None

    m = folium.Map(location=[20, 0], zoom_start=2,
                   tiles="CartoDB positron")

    # Aggregate sequences per country
    country_counts = (
        geo_df.groupby("country")["sequences_count"]
              .sum()
              .reset_index()
    )

    # First detection markers per variant per country
    first_detections = (
        geo_df.sort_values("first_detection")
              .drop_duplicates(subset=["who_label"])
    )

    # Country centroids (subset)
    centroids = {
        "United Kingdom": (51.5, -0.1),
        "India":          (20.6, 78.9),
        "South Africa":   (-30.6, 22.9),
        "United States":  (37.1, -95.7),
        "Brazil":         (-14.2, -51.9),
        "Denmark":        (56.3, 9.5),
        "Germany":        (51.2, 10.5),
        "France":         (46.2, 2.2),
        "China":          (35.9, 104.2),
        "Australia":      (-25.3, 133.8),
        "Japan":          (36.2, 138.3),
        "Singapore":      (1.4, 103.8),
        "Kenya":          (-0.0, 37.9),
        "Nigeria":        (9.1, 8.7),
    }

    mc = MarkerCluster().add_to(m)

    for _, row in first_detections.iterrows():
        country = row.get("country", "Unknown")
        if country not in centroids:
            continue
        lat, lon = centroids[country]
        variant  = row["who_label"]
        colour   = WHO_COLOURS.get(variant, "#888")
        date_str = str(row["first_detection"])[:10]

        folium.CircleMarker(
            location=[lat, lon],
            radius=12,
            color=colour,
            fill=True,
            fill_color=colour,
            fill_opacity=0.8,
            popup=folium.Popup(
                f"<b>{variant}</b><br>"
                f"First detected: {date_str}<br>"
                f"Country: {country}",
                max_width=200,
            ),
            tooltip=f"{variant} — {country} ({date_str})",
        ).add_to(mc)

    return m


def geographic_plotly(geo_df: pd.DataFrame) -> go.Figure:
    """
    Plotly choropleth (works without folium).
    Shows log-scaled sequence counts per country.
    """
    country_counts = (
        geo_df.groupby("country")["sequences_count"]
              .sum()
              .reset_index()
    )

    fig = px.choropleth(
        country_counts,
        locations="country",
        locationmode="country names",
        color="sequences_count",
        color_continuous_scale="Blues",
        title="SARS-CoV-2 sequences per country",
        labels={"sequences_count": "Sequences"},
    )

    # Overlay first-detection scatter
    variant_origins = (
        geo_df.sort_values("first_detection")
              .drop_duplicates(subset=["who_label"])
    )

    country_coords = {
        "United Kingdom": (51.5, -0.1),
        "India":          (20.6, 78.9),
        "South Africa":   (-30.6, 22.9),
        "United States":  (37.1, -95.7),
        "Brazil":         (-14.2, -51.9),
        "Denmark":        (56.3, 9.5),
    }

    lats, lons, texts, colours = [], [], [], []
    for _, row in variant_origins.iterrows():
        c = row.get("country", "")
        if c in country_coords:
            lat, lon = country_coords[c]
            lats.append(lat)
            lons.append(lon)
            texts.append(f"{row['who_label']} — first detected {str(row['first_detection'])[:7]}")
            colours.append(WHO_COLOURS.get(row["who_label"], "#888"))

    fig.add_trace(go.Scattergeo(
        lat=lats, lon=lons,
        mode="markers+text",
        marker=dict(size=14, color=colours, symbol="star"),
        text=[t.split("—")[0].strip() for t in texts],
        textposition="top center",
        hovertext=texts,
        hoverinfo="text",
        name="Variant origin",
    ))

    _apply_theme(fig)
    return fig


# ── 3. Spike protein linear map ────────────────────────────────────────────

def spike_protein_map(
    hotspots: pd.DataFrame,
    highlight_variant: Optional[str] = None,
) -> go.Figure:
    """
    Horizontal linear map of the 1,273-aa spike protein with:
      • Coloured domain blocks
      • Hotspot circles scaled by mutation count
      • Key mutation labels
    """
    fig = go.Figure()
    SPIKE_LEN = 1273

    # Draw domain rectangles
    for domain, (start, end) in SPIKE_DOMAINS.items():
        col = DOMAIN_COLOURS.get(domain, "#ccc")
        fig.add_shape(
            type="rect",
            x0=start, x1=end, y0=0.3, y1=0.7,
            fillcolor=col, opacity=0.85,
            line=dict(color="white", width=1),
        )
        mid = (start + end) / 2
        if (end - start) > 30:
            fig.add_annotation(
                x=mid, y=0.5,
                text=domain if (end-start) > 60 else domain[:3],
                showarrow=False,
                font=dict(size=9, color="#333"),
            )

    # Plot hotspot circles
    if not hotspots.empty:
        hs = hotspots.copy()
        hs["size"] = np.sqrt(hs["mutation_count"]) * 3 + 4
        hs["size"] = hs["size"].clip(upper=20)

        domain_col = hs["domain"].map(DOMAIN_COLOURS).fillna("#888")

        fig.add_trace(go.Scatter(
            x=hs["position"],
            y=[0.85] * len(hs),
            mode="markers",
            marker=dict(
                size=hs["size"],
                color=domain_col,
                opacity=0.8,
                line=dict(color="white", width=1),
            ),
            text=hs["top_mutations"],
            hovertemplate=(
                "<b>Position %{x}</b><br>"
                "Mutations: %{text}<br>"
                "Count: %{customdata}<extra></extra>"
            ),
            customdata=hs["mutation_count"],
            name="Hotspots",
        ))

    # Key mutation labels
    key_pos = {m[0]: m[1] for m in KEY_MUTATIONS}
    for mut, pos in key_pos.items():
        fig.add_annotation(
            x=pos, y=1.05,
            text=mut,
            showarrow=True,
            arrowhead=2,
            arrowsize=0.8,
            arrowwidth=1,
            arrowcolor="#555",
            font=dict(size=8, color="#333"),
            ax=0, ay=-15,
        )

    fig.update_layout(
        title="SARS-CoV-2 spike protein — mutational hotspots",
        xaxis=dict(
            range=[0, SPIKE_LEN],
            title="Amino acid position",
            showgrid=False,
        ),
        yaxis=dict(range=[0, 1.2], showticklabels=False, showgrid=False),
        height=250,
        showlegend=False,
    )
    _apply_theme(fig)
    return fig


def spike_rbd_detail(hotspots: pd.DataFrame) -> go.Figure:
    """Zoomed-in bar chart of RBD positions (306-527) mutation counts."""
    rbd = hotspots[hotspots["domain"].isin(["RBD", "RBM"])].copy()
    if rbd.empty:
        return _empty_fig("No RBD hotspot data")

    rbd = rbd.sort_values("position")
    colours = rbd["domain"].map({"RBD": "#F09595", "RBM": "#E24B4A"}).fillna("#ccc")

    fig = go.Figure(go.Bar(
        x=rbd["position"],
        y=rbd["mutation_count"],
        text=rbd["top_mutations"],
        marker_color=colours,
        hovertemplate="<b>Position %{x}</b><br>Count: %{y}<br>%{text}<extra></extra>",
    ))
    fig.update_layout(
        title="RBD mutation frequency by position (residues 306–527)",
        xaxis_title="Amino acid position",
        yaxis_title="Mutation count",
        xaxis=dict(range=[306, 527]),
    )
    _apply_theme(fig)
    return fig


# ── 4. Phylogenetic tree ───────────────────────────────────────────────────

def phylogenetic_tree_plotly() -> go.Figure:
    """
    Build a Plotly-based phylogenetic tree from the known variant hierarchy.
    Layout is computed manually; for large trees, use ETE3 (see ete3_tree()).
    """
    nodes = [
        # id, x, y, who_label, defining_muts
        ("Wuhan-Hu-1",     0,  0.00, "Ancestral",  "Reference genome"),
        ("B.1 clade",      1,  0.00, "Ancestral",  "D614G backbone"),
        ("Alpha B.1.1.7",  2,  0.40, "Alpha",      "N501Y, P681H"),
        ("Beta B.1.351",   2,  0.20, "Beta",       "K417N, E484K, N501Y"),
        ("Gamma P.1",      2,  0.00, "Gamma",      "K417T, E484K, N501Y"),
        ("Delta B.1.617.2",2, -0.20, "Delta",      "L452R, T478K, P681R"),
        ("AY.4 Delta+",    3, -0.25, "Delta",      "+K417N"),
        ("Omicron BA.1",   3,  0.55, "Omicron",    "30+ spike muts"),
        ("BA.2",           4,  0.65, "Omicron",    "28 spike muts"),
        ("BA.4/BA.5",      4,  0.45, "Omicron",    "L452R, F486V"),
        ("XBB",            5,  0.75, "Omicron",    "Recombinant"),
        ("XBB.1.5",        6,  0.80, "Omicron",    "F486P, N460K"),
        ("EG.5 / HV.1",   7,  0.82, "Omicron",    "F456L"),
        ("BA.2.86",        5,  0.60, "Omicron",    "+30 muts from BA.2"),
        ("JN.1",           6,  0.55, "Omicron",    "L455S"),
    ]

    edges = [
        ("Wuhan-Hu-1",      "B.1 clade"),
        ("B.1 clade",       "Alpha B.1.1.7"),
        ("B.1 clade",       "Beta B.1.351"),
        ("B.1 clade",       "Gamma P.1"),
        ("B.1 clade",       "Delta B.1.617.2"),
        ("Delta B.1.617.2", "AY.4 Delta+"),
        ("B.1 clade",       "Omicron BA.1"),
        ("Omicron BA.1",    "BA.2"),
        ("Omicron BA.1",    "BA.4/BA.5"),
        ("BA.2",            "XBB"),
        ("XBB",             "XBB.1.5"),
        ("XBB.1.5",         "EG.5 / HV.1"),
        ("BA.2",            "BA.2.86"),
        ("BA.2.86",         "JN.1"),
    ]

    node_dict = {n[0]: n for n in nodes}

    # Edge traces
    edge_x, edge_y = [], []
    for src, dst in edges:
        sx, sy = node_dict[src][1], node_dict[src][2]
        dx, dy = node_dict[dst][1], node_dict[dst][2]
        edge_x += [sx, dx, None]
        edge_y += [sy, dy, None]

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=edge_x, y=edge_y,
        mode="lines",
        line=dict(width=1.5, color="#B4B2A9"),
        hoverinfo="none",
        showlegend=False,
    ))

    # Node traces by WHO label
    who_groups: Dict[str, list] = {}
    for n in nodes:
        who_groups.setdefault(n[3], []).append(n)

    for who, group in who_groups.items():
        xs = [g[1] for g in group]
        ys = [g[2] for g in group]
        names  = [g[0] for g in group]
        muts   = [g[4] for g in group]

        fig.add_trace(go.Scatter(
            x=xs, y=ys,
            mode="markers+text",
            name=who,
            text=names,
            textposition="middle right",
            marker=dict(size=14, color=WHO_COLOURS.get(who, "#888"),
                        line=dict(color="white", width=2)),
            customdata=muts,
            hovertemplate="<b>%{text}</b><br>Defining muts: %{customdata}<extra></extra>",
        ))

    fig.update_layout(
        title="SARS-CoV-2 variant phylogenetic tree (simplified)",
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False,
                   title="Divergence (schematic)"),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        height=500,
        legend=dict(orientation="v", x=0.02, y=0.98),
    )
    _apply_theme(fig)
    return fig


# ── 5. Co-occurrence heatmap ───────────────────────────────────────────────

def cooccurrence_heatmap(cooccur_df: pd.DataFrame, top_n: int = 15) -> go.Figure:
    """Symmetric heatmap of mutation co-occurrence odds ratios."""
    if cooccur_df.empty:
        return _empty_fig("No co-occurrence data")

    top = cooccur_df.head(top_n * 5)
    muts = list(set(top["mut_a"].tolist() + top["mut_b"].tolist()))[:top_n]

    matrix = pd.DataFrame(1.0, index=muts, columns=muts)
    for _, row in top.iterrows():
        if row["mut_a"] in muts and row["mut_b"] in muts:
            matrix.loc[row["mut_a"], row["mut_b"]] = row["odds_ratio"]
            matrix.loc[row["mut_b"], row["mut_a"]] = row["odds_ratio"]

    fig = go.Figure(go.Heatmap(
        z=matrix.values,
        x=muts, y=muts,
        colorscale="RdBu",
        zmid=1.0,
        colorbar=dict(title="Odds ratio"),
        hovertemplate="<b>%{y} × %{x}</b><br>Odds ratio: %{z:.2f}<extra></extra>",
    ))
    fig.update_layout(
        title="Mutation co-occurrence odds ratios (spike protein)",
        height=450,
    )
    _apply_theme(fig)
    return fig


# ── 6. Sweep speed bar chart ───────────────────────────────────────────────

def sweep_speed_chart(sweep_df: pd.DataFrame) -> go.Figure:
    """Horizontal bar chart showing days to sweep from 5% to 50% frequency."""
    if sweep_df.empty or "sweep_days" not in sweep_df.columns:
        return _empty_fig("No sweep data")

    df = sweep_df.dropna(subset=["sweep_days"]).sort_values("sweep_days")
    colours = [FUNCTIONAL_COLOURS.get(
        next((m[2] for m in KEY_MUTATIONS if m[0] == mut), "Other"), "#888")
        for mut in df["mutation"]]

    fig = go.Figure(go.Bar(
        y=df["mutation"],
        x=df["sweep_days"],
        orientation="h",
        marker_color=colours,
        text=df["sweep_days"].astype(int).astype(str) + " days",
        textposition="outside",
        hovertemplate="<b>%{y}</b><br>Sweep time: %{x} days<extra></extra>",
    ))
    fig.update_layout(
        title="Mutation sweep speed (5% → 50% global frequency)",
        xaxis_title="Days",
        yaxis_title="",
        height=max(300, len(df) * 30 + 80),
    )
    _apply_theme(fig)
    return fig


# ── Shared theme ───────────────────────────────────────────────────────────

def _apply_theme(fig: go.Figure):
    fig.update_layout(
        font=dict(family="Arial, sans-serif", size=12),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        margin=dict(l=50, r=20, t=50, b=50),
        xaxis=dict(gridcolor="rgba(128,128,128,0.15)"),
        yaxis=dict(gridcolor="rgba(128,128,128,0.15)"),
    )


def _empty_fig(msg: str) -> go.Figure:
    fig = go.Figure()
    fig.add_annotation(text=msg, x=0.5, y=0.5, showarrow=False,
                       xref="paper", yref="paper", font=dict(size=14))
    _apply_theme(fig)
    return fig
