"""
app.py — COVID-19 Mutation Tracker
───────────────────────────────────
Main Streamlit dashboard.

Run:
    streamlit run app.py

First run will auto-generate synthetic sample data.
To use real data, run the pipeline first:
    python pipeline/fetch_data.py --source nextstrain --max-sequences 10000
    python pipeline/alignment.py
"""

import json
import sys
from pathlib import Path

import pandas as pd
import streamlit as st

# ── Page config (must be first Streamlit call) ─────────────────────────────
st.set_page_config(
    page_title="SARS-CoV-2 Mutation Tracker",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

sys.path.insert(0, str(Path(__file__).parent))
from config import OUTPUT_DIR, KEY_MUTATIONS, VARIANT_DEFINITIONS
from visualization.plots import (
    timeline_plot,
    geographic_plotly,
    spike_protein_map,
    spike_rbd_detail,
    sweep_speed_chart,
    WHO_COLOURS,
)


# ── Cached data loaders ────────────────────────────────────────────────────

@st.cache_data(show_spinner=False)
def load_all_data():
    """Load or generate all processed data files."""
    mut_table_path = OUTPUT_DIR / "mutation_table.parquet"

    # Auto-run pipeline if data not present
    if not mut_table_path.exists():
        with st.spinner("🔬 Generating sample data (first run — takes ~10s) …"):
            from pipeline.fetch_data   import _generate_sample_metadata
            from pipeline.alignment    import build_mutation_table
            from analysis.frequency_analysis import run_all_analyses

            meta_path = _generate_sample_metadata(n=5000)
            df        = build_mutation_table(meta_path)
            run_all_analyses(df)

    data = {}

    def _load(name):
        p = OUTPUT_DIR / f"{name}.parquet"
        return pd.read_parquet(p) if p.exists() else pd.DataFrame()

    data["mutations"]        = _load("mutation_table")
    data["variant_freq"]     = _load("variant_frequency")
    data["mutation_freq"]    = _load("mutation_frequency_time")
    data["convergent"]       = _load("convergent_evolution")
    data["hotspots"]         = _load("hotspots")
    data["sweeps"]           = _load("sweep_speed")
    data["geo"]              = _load("geographic_spread")

    # Fix period columns
    for key in ["variant_freq", "mutation_freq"]:
        df = data[key]
        if not df.empty and "year_month" in df.columns:
            if not hasattr(df["year_month"].dtype, "freq"):
                df["year_month"] = pd.PeriodIndex(df["year_month"], freq="M")
            data[key] = df

    return data


# ── Sidebar ────────────────────────────────────────────────────────────────

def render_sidebar(data):
    st.sidebar.image(
        "https://upload.wikimedia.org/wikipedia/commons/thumb/8/82/"
        "SARS-CoV-2_without_background.png/220px-SARS-CoV-2_without_background.png",
        width=140,
    )
    st.sidebar.title("🧬 Filters")

    # Date range
    mutations_df = data["mutations"]
    if not mutations_df.empty and "date" in mutations_df.columns:
        mutations_df["date"] = pd.to_datetime(mutations_df["date"])
        min_date = mutations_df["date"].min().date()
        max_date = mutations_df["date"].max().date()
        date_range = st.sidebar.date_input(
            "Date range",
            value=(min_date, max_date),
            min_value=min_date,
            max_value=max_date,
        )
    else:
        date_range = None

    # Region filter
    if not mutations_df.empty and "region" in mutations_df.columns:
        regions = ["All"] + sorted(mutations_df["region"].dropna().unique().tolist())
        region = st.sidebar.selectbox("Region", regions)
    else:
        region = "All"

    # Variant filter
    if not mutations_df.empty and "who_label" in mutations_df.columns:
        variants = ["All"] + sorted(mutations_df["who_label"].dropna().unique().tolist())
        variant = st.sidebar.selectbox("Variant (WHO label)", variants)
    else:
        variant = "All"

    st.sidebar.divider()
    st.sidebar.markdown("**Data source**")
    source_label = "Synthetic sample data (5,000 sequences)"
    real_path = OUTPUT_DIR / "nextstrain_metadata.parquet"
    if real_path.exists():
        df_size = len(pd.read_parquet(real_path))
        if df_size > 5000:
            source_label = f"Nextstrain open data ({df_size:,} sequences)"
    st.sidebar.info(source_label)

    return date_range, region, variant


# ── Apply filters ─────────────────────────────────────────────────────────

def apply_filters(data, date_range, region, variant):
    df = data["mutations"].copy()
    if df.empty:
        return data

    df["date"] = pd.to_datetime(df["date"])

    if date_range and len(date_range) == 2:
        df = df[
            (df["date"] >= pd.Timestamp(date_range[0])) &
            (df["date"] <= pd.Timestamp(date_range[1]))
        ]
    if region != "All" and "region" in df.columns:
        df = df[df["region"] == region]
    if variant != "All" and "who_label" in df.columns:
        df = df[df["who_label"] == variant]

    filtered = dict(data)
    filtered["mutations"] = df
    return filtered


# ── Metric cards ───────────────────────────────────────────────────────────

def render_metrics(df):
    c1, c2, c3, c4, c5 = st.columns(5)
    with c1:
        st.metric("Sequences", f"{len(df):,}")
    with c2:
        n_muts = sum(
            len(json.loads(m)) for m in df["spike_mutations"].dropna()
            if isinstance(m, str)
        )
        unique = set()
        for m in df["spike_mutations"].dropna():
            if isinstance(m, str):
                unique.update(json.loads(m))
        st.metric("Unique mutations", f"{len(unique):,}")
    with c3:
        n_lineages = df["who_label"].nunique() if "who_label" in df.columns else 0
        st.metric("WHO variants", n_lineages)
    with c4:
        n_countries = df["country"].nunique() if "country" in df.columns else 0
        st.metric("Countries", n_countries)
    with c5:
        avg_rbd = df["rbd_mut_count"].mean() if "rbd_mut_count" in df.columns else 0
        st.metric("Avg RBD mutations", f"{avg_rbd:.1f}")


# ── Tab: Overview ─────────────────────────────────────────────────────────

def tab_overview(data, filtered):
    st.subheader("Genomic surveillance overview")
    render_metrics(filtered["mutations"])

    st.divider()

    col_left, col_right = st.columns([3, 2])

    with col_left:
        vf = data["variant_freq"]
        if not vf.empty:
            fig = timeline_plot(vf, mode="variant")
            st.plotly_chart(fig, width='stretch', key="chart_1")

    with col_right:
        st.markdown("**Key evolutionary events**")
        events = [
            ("2020-02", "D614G",  "Transmissibility", "D614G sweeps globally — now universal backbone"),
            ("2020-09", "Alpha",  "Alpha",  "N501Y + P681H confer ~50% transmission advantage"),
            ("2020-10", "Beta",   "Beta",   "K417N + E484K: immune escape first documented"),
            ("2020-10", "E484K",  "Immune", "Convergent evolution: E484K in 3 independent lineages"),
            ("2021-04", "Delta",  "Delta",  "L452R + T478K; P681R strengthens furin cleavage"),
            ("2021-11", "Omicron","Omicron","30+ spike mutations; unprecedented immune evasion"),
            ("2022-10", "XBB",   "Omicron","Recombinant lineage; F486P key for ACE2 + escape"),
            ("2023-08", "JN.1",  "Omicron","BA.2.86 descendant; L455S novel ACE2 binding shift"),
        ]
        for date, label, variant, desc in events:
            col = WHO_COLOURS.get(variant, "#888")
            st.markdown(
                f"<div style='border-left:3px solid {col}; "
                f"padding:6px 10px; margin:4px 0; border-radius:0 4px 4px 0'>"
                f"<span style='font-size:11px; color:#888'>{date}</span> "
                f"<strong>{label}</strong><br>"
                f"<span style='font-size:12px'>{desc}</span></div>",
                unsafe_allow_html=True,
            )




# ── Tab: Mutation Timeline ─────────────────────────────────────────────────

def tab_timeline(data):
    st.subheader("Mutation frequency over time")

    mode = st.radio("View mode", ["Variant dominance", "Individual mutations"],
                    horizontal=True)

    if mode == "Individual mutations":
        all_muts = [m[0] for m in KEY_MUTATIONS]
        selected = st.multiselect(
            "Select mutations",
            all_muts,
            default=["D614G", "N501Y", "E484K", "L452R", "P681H"],
        )
        mf = data["mutation_freq"]
        fig = timeline_plot(data["variant_freq"], mutation_freq=mf,
                            mode="mutation", selected_mutations=selected)
    else:
        fig = timeline_plot(data["variant_freq"], mode="variant")

    st.plotly_chart(fig, width='stretch', key="chart_2")

    # Sweep speed
    st.subheader("Selective sweep speed")
    st.caption("Days for each mutation to rise from 5% → 50% global frequency")
    sweeps = data["sweeps"]
    if not sweeps.empty:
        col1, col2 = st.columns([2, 1])
        with col1:
            st.plotly_chart(sweep_speed_chart(sweeps), key="chart_3", width='stretch')
        with col2:
            st.dataframe(
                sweeps[["mutation", "first_seen", "reached_50pct", "sweep_days", "peak_freq"]]
                      .rename(columns={
                          "mutation":     "Mutation",
                          "first_seen":   "First seen",
                          "reached_50pct":"Reached 50%",
                          "sweep_days":   "Sweep (days)",
                          "peak_freq":    "Peak freq %",
                      }),
                width='stretch',
                hide_index=True,
            )


# ── Tab: Geographic ───────────────────────────────────────────────────────

def tab_geographic(data):
    st.subheader("Geographic spread & variant origins")

    geo = data["geo"]
    if geo.empty:
        st.warning("No geographic data available.")
        return

    # Choropleth
    fig = geographic_plotly(geo)
    st.plotly_chart(fig, width='stretch', key="chart_4")

    col1, col2 = st.columns(2)

    with col1:
        st.markdown("**First detection by variant**")
        first = (
            geo.sort_values("first_detection")
               .drop_duplicates(subset=["who_label"])
               [["who_label", "country", "first_detection", "sequences_count"]]
               .rename(columns={
                   "who_label":       "Variant",
                   "country":         "First country",
                   "first_detection": "First detected",
                   "sequences_count": "Sequences",
               })
        )
        st.dataframe(first, width='stretch', hide_index=True)

    with col2:
        st.markdown("**Top contributing countries**")
        top_countries = (
            geo.groupby("country")["sequences_count"]
               .sum()
               .reset_index()
               .sort_values("sequences_count", ascending=False)
               .head(15)
               .rename(columns={"country": "Country", "sequences_count": "Sequences"})
        )
        st.dataframe(top_countries, width='stretch', hide_index=True)


# ── Tab: Spike Protein ────────────────────────────────────────────────────

def tab_spike(data):
    st.subheader("Spike protein mutational landscape")

    hs = data["hotspots"]

    # Full protein map
    st.markdown("**Full spike protein (1,273 aa) — hotspot map**")
    st.caption("Circle size = number of sequences with a mutation at that position")
    fig_full = spike_protein_map(hs)
    st.plotly_chart(fig_full, width='stretch', key="chart_5")

    col1, col2 = st.columns(2)

    with col1:
        st.markdown("**RBD zoom (positions 306–527)**")
        fig_rbd = spike_rbd_detail(hs)
        st.plotly_chart(fig_rbd, width='stretch', key="chart_6")

    with col2:
        st.markdown("**Key mutation functional categories**")
        import plotly.graph_objects as go
        from config import FUNCTIONAL_COLOURS
        cats = {}
        for _, _, cat, _ in KEY_MUTATIONS:
            cats[cat] = cats.get(cat, 0) + 1
        fig_pie = go.Figure(go.Pie(
            labels=list(cats.keys()),
            values=list(cats.values()),
            marker_colors=[FUNCTIONAL_COLOURS.get(c, "#888") for c in cats],
            hole=0.4,
            textinfo="label+percent",
        ))
        fig_pie.update_layout(
            showlegend=False,
            margin=dict(l=10, r=10, t=30, b=10),
            paper_bgcolor="rgba(0,0,0,0)",
        )
        st.plotly_chart(fig_pie, width='stretch', key="chart_7")

    # Domain breakdown table
    st.markdown("**Mutation count per domain**")
    if not hs.empty:
        domain_summary = (
            hs.groupby("domain")
              .agg(
                  total_mutations=("mutation_count", "sum"),
                  positions=("position", "count"),
                  top_positions=("position", lambda x: ", ".join(x.nlargest(3).astype(str)))
              )
              .sort_values("total_mutations", ascending=False)
              .reset_index()
              .rename(columns={
                  "domain":           "Domain",
                  "total_mutations":  "Total mutations",
                  "positions":        "Mutated positions",
                  "top_positions":    "Top positions",
              })
        )
        st.dataframe(domain_summary, width='stretch', hide_index=True)


# ── Tab: Convergent Evolution ─────────────────────────────────────────────

def tab_convergent(data):
    st.subheader("Convergent evolution")
    st.caption("Mutations that appeared independently in multiple distinct variants — a sign of strong selective pressure.")

    conv = data["convergent"]
    if conv.empty:
        st.info("No convergent evolution data available.")
        return

    import plotly.express as px
    top_conv = conv.nlargest(12, "appearance_count")
    fig = px.bar(
        top_conv,
        x="mutation",
        y="appearance_count",
        color="appearance_count",
        color_continuous_scale="Reds",
        labels={"mutation": "Mutation", "appearance_count": "Number of variants"},
        title="Spike mutations with strongest convergent evolution",
    )
    fig.update_layout(
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        showlegend=False,
        coloraxis_showscale=False,
        xaxis_title="Mutation",
        yaxis_title="Independent variant appearances",
    )
    st.plotly_chart(fig, width='stretch', key="chart_8")

    st.markdown("**What this means:** Each bar shows a spike mutation that evolved independently in multiple SARS-CoV-2 lineages. When the same mutation arises repeatedly across unrelated variants, it indicates the virus is under strong evolutionary pressure to acquire that change — usually for increased transmissibility or immune escape.")

    col1, col2 = st.columns(2)
    with col1:
        st.markdown("**Top convergent mutations**")
        st.dataframe(
            conv.nlargest(8, "appearance_count")[["mutation", "appearance_count", "lineages"]]
                .rename(columns={
                    "mutation":         "Mutation",
                    "appearance_count": "Variants",
                    "lineages":         "Found in",
                }),
            width='stretch',
            hide_index=True,
        )
    with col2:
        st.markdown("**Key example: E484K**")
        st.info(
            "E484K is the classic example of convergent evolution in SARS-CoV-2. "
            "It appeared independently in the Beta (South Africa), Gamma (Brazil), "
            "and Iota (USA) lineages — all acquiring the same immune escape mutation "
            "without sharing a common ancestor."
        )





# ── Tab: Data Explorer ─────────────────────────────────────────────────────

def tab_data(filtered):
    st.subheader("Data explorer")

    df = filtered["mutations"].copy()
    if df.empty:
        st.warning("No data matches current filters.")
        return

    # Search
    search = st.text_input("Search by strain ID, country, or mutation",
                            placeholder="e.g. EPI_ISL_402125 or N501Y …")
    if search:
        mask = (
            df["strain"].str.contains(search, case=False, na=False) |
            df["country"].str.contains(search, case=False, na=False) |
            df["spike_mutations"].str.contains(search, case=False, na=False)
        )
        df = df[mask]

    # Column selection
    show_cols = ["strain", "date", "country", "region", "lineage", "who_label",
                 "spike_mut_count", "rbd_mut_count"]
    show_cols = [c for c in show_cols if c in df.columns]

    st.dataframe(
        df[show_cols]
          .rename(columns={
              "strain":          "Sample ID",
              "date":            "Date",
              "country":         "Country",
              "region":          "Region",
              "lineage":         "Pango lineage",
              "who_label":       "WHO label",
              "spike_mut_count": "Spike mutations",
              "rbd_mut_count":   "RBD mutations",
          })
          .sort_values("Date", ascending=False)
          .head(500),
        width='stretch',
        hide_index=True,
    )

    st.caption(f"Showing up to 500 of {len(df):,} filtered sequences")

    # Download
    csv_data = df[show_cols].to_csv(index=False)
    st.download_button(
        label="⬇ Download filtered data (CSV)",
        data=csv_data,
        file_name="covid_mutations_filtered.csv",
        mime="text/csv",
    )


# ── Tab: About ────────────────────────────────────────────────────────────

def tab_about():
    st.subheader("About")
    st.markdown("""
## SARS-CoV-2 Mutation Tracker

A full genomic-surveillance pipeline that tracks how SARS-CoV-2 spike protein
mutations emerge, spread geographically, and affect viral fitness over time.

### Pipeline
| Module | Script | What it does |
|--------|--------|-------------|
| Data fetch | `pipeline/fetch_data.py` | Downloads NCBI reference + Nextstrain/GISAID metadata |
| Alignment | `pipeline/alignment.py` | Parses/aligns sequences, calls amino-acid mutations |
| Analysis | `analysis/frequency_analysis.py` | Frequency, co-occurrence, convergence, hotspots, sweeps |
| Visualisation | `visualization/plots.py` | All Plotly/Folium figures |
| Dashboard | `app.py` | This Streamlit app |

### Data sources
- **Nextstrain open** — pre-processed metadata with mutation annotations (no account required)
- **GISAID** — full sequences + metadata (free account required at [gisaid.org](https://gisaid.org))
- **NCBI** — Wuhan-Hu-1 reference genome MN908947

### Key analyses
1. **Variant dominance curves** — monthly WHO-labelled variant frequency
2. **Mutation sweep speed** — days from 5% → 50% global frequency
3. **Convergent evolution** — mutations appearing independently in ≥2 lineages
4. **Co-occurrence / LD** — chi-squared test + odds ratio for mutation linkage
5. **Spike hotspot map** — positional mutation density across all 1,273 aa
6. **Geographic spread** — first-detection date per variant per country

### References
- [Nextstrain SARS-CoV-2](https://nextstrain.org/ncov/open/global)
- [GISAID Initiative](https://gisaid.org)
- [CoV-Spectrum](https://cov-spectrum.org)
- [Outbreak.info](https://outbreak.info)

### Citation
> This dashboard uses publicly available SARS-CoV-2 genomic data.
> We gratefully acknowledge all data contributors.
    """)


# ── Main ───────────────────────────────────────────────────────────────────

def main():
    # Header
    st.markdown(
        "<h1 style='margin-bottom:0'>🧬 SARS-CoV-2 Mutation Tracker</h1>"
        "<p style='color:#888; margin-top:4px'>Tracking how SARS-CoV-2 spike mutations emerge, spread, and affect viral fitness · Nextstrain open data</p>",
        unsafe_allow_html=True,
    )
    st.divider()

    # Load data
    with st.spinner("Loading data …"):
        data = load_all_data()

    # Sidebar filters
    date_range, region, variant = render_sidebar(data)

    # Apply filters
    filtered = apply_filters(data, date_range, region, variant)

    # Tabs
    tabs = st.tabs([
        "📊 Overview",
        "📈 Mutation Timeline",
        "🌍 Geographic Spread",
        "🔬 Spike Protein",
        "🧬 Convergent Evolution",
        "🗂 Data Explorer",
        "ℹ About",
    ])

    with tabs[0]:
        tab_overview(data, filtered)
    with tabs[1]:
        tab_timeline(data)
    with tabs[2]:
        tab_geographic(data)
    with tabs[3]:
        tab_spike(data)
    with tabs[4]:
        tab_convergent(data)
    with tabs[5]:
        tab_data(filtered)
    with tabs[6]:
        tab_about()


if __name__ == "__main__":
    main()
