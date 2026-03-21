"""
app.py — SARS-CoV-2 Mutation Tracker
"""

import json
import sys
from pathlib import Path

import pandas as pd
import streamlit as st

st.set_page_config(
    page_title="SARS-CoV-2 // Mutation Tracker",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@400;600&family=Inter:wght@300;400;500;600&display=swap');
html, body, [class*="css"] { font-family: 'Inter', sans-serif; }
.stApp { background-color: #0d1117; color: #e6edf3; }
[data-testid="stSidebar"] { background-color: #090c10 !important; border-right: 1px solid #21262d; }
[data-testid="stSidebar"] * { color: #8b949e !important; }
.stTabs [data-baseweb="tab-list"] { background-color: #0d1117; border-bottom: 1px solid #21262d; gap: 0px; }
.stTabs [data-baseweb="tab"] { background-color: transparent; color: #8b949e; font-size: 12px; font-family: 'JetBrains Mono', monospace; padding: 8px 16px; border: none; border-bottom: 2px solid transparent; }
.stTabs [aria-selected="true"] { color: #58a6ff !important; border-bottom: 2px solid #58a6ff !important; background-color: transparent !important; }
.stTabs [data-baseweb="tab"]:hover { color: #e6edf3 !important; background-color: #161b22 !important; }
[data-testid="stMetric"] { background-color: #161b22; border: 1px solid #21262d; border-radius: 6px; padding: 16px 20px; }
[data-testid="stMetricLabel"] { font-size: 11px !important; text-transform: uppercase; letter-spacing: 0.08em; color: #8b949e !important; font-family: 'JetBrains Mono', monospace !important; }
[data-testid="stMetricValue"] { font-size: 26px !important; font-weight: 600 !important; color: #e6edf3 !important; font-family: 'JetBrains Mono', monospace !important; }
[data-testid="stDataFrame"] { border: 1px solid #21262d !important; border-radius: 6px; }
.stButton > button, .stDownloadButton > button { background-color: #21262d !important; color: #58a6ff !important; border: 1px solid #30363d !important; border-radius: 6px !important; font-family: 'JetBrains Mono', monospace !important; font-size: 12px !important; }
.stSelectbox > div > div, .stTextInput > div > div > input { background-color: #161b22 !important; border: 1px solid #30363d !important; color: #e6edf3 !important; border-radius: 6px !important; }
hr { border-color: #21262d !important; margin: 16px 0; }
.stAlert { background-color: #161b22 !important; border: 1px solid #21262d !important; border-radius: 6px !important; color: #8b949e !important; }
h1, h2, h3 { color: #e6edf3 !important; font-weight: 600 !important; }
::-webkit-scrollbar { width: 6px; height: 6px; }
::-webkit-scrollbar-track { background: #0d1117; }
::-webkit-scrollbar-thumb { background: #30363d; border-radius: 3px; }
</style>
""", unsafe_allow_html=True)

sys.path.insert(0, str(Path(__file__).parent))
from config import OUTPUT_DIR, KEY_MUTATIONS, VARIANT_DEFINITIONS
from visualization.plots import (
    timeline_plot, geographic_plotly, spike_protein_map,
    spike_rbd_detail, sweep_speed_chart, WHO_COLOURS,
)

PLOTLY_DARK = dict(
    paper_bgcolor="#0d1117", plot_bgcolor="#0d1117",
    font=dict(color="#8b949e", family="JetBrains Mono, monospace", size=11),
    xaxis=dict(gridcolor="#21262d", linecolor="#30363d"),
    yaxis=dict(gridcolor="#21262d", linecolor="#30363d"),
    legend=dict(bgcolor="#161b22", bordercolor="#30363d", borderwidth=1),
)

def dark_fig(fig):
    fig.update_layout(**PLOTLY_DARK)
    return fig

def section(title, subtitle=""):
    sub = f"<br><span style='font-size:12px; color:#6e7681'>{subtitle}</span>" if subtitle else ""
    st.markdown(
        f"<div style='margin:8px 0 14px'><span style='font-family:JetBrains Mono,monospace;"
        f"font-size:11px;color:#6e7681;text-transform:uppercase;letter-spacing:0.1em'>"
        f"{title}</span>{sub}</div>", unsafe_allow_html=True)

@st.cache_data(show_spinner=False)
def load_all_data():
    mut_table_path = OUTPUT_DIR / "mutation_table.parquet"
    if not mut_table_path.exists():
        with st.spinner("initialising..."):
            from pipeline.fetch_data import _generate_sample_metadata
            from pipeline.alignment import build_mutation_table
            from analysis.frequency_analysis import run_all_analyses
            meta_path = _generate_sample_metadata(n=5000)
            df = build_mutation_table(meta_path)
            run_all_analyses(df)
    data = {}
    def _load(name):
        p = OUTPUT_DIR / f"{name}.parquet"
        return pd.read_parquet(p) if p.exists() else pd.DataFrame()
    data["mutations"]     = _load("mutation_table")
    data["variant_freq"]  = _load("variant_frequency")
    data["mutation_freq"] = _load("mutation_frequency_time")
    data["convergent"]    = _load("convergent_evolution")
    data["hotspots"]      = _load("hotspots")
    data["sweeps"]        = _load("sweep_speed")
    data["geo"]           = _load("geographic_spread")
    for key in ["variant_freq", "mutation_freq"]:
        df = data[key]
        if not df.empty and "year_month" in df.columns:
            if not hasattr(df["year_month"].dtype, "freq"):
                df["year_month"] = pd.PeriodIndex(df["year_month"], freq="M")
            data[key] = df
    return data

def render_sidebar(data):
    st.sidebar.markdown(
        "<div style='padding:16px 0 8px;font-family:JetBrains Mono,monospace;"
        "font-size:13px;color:#58a6ff;letter-spacing:0.05em'>// SARS-CoV-2<br>"
        "<span style='font-size:11px;color:#6e7681'>mutation tracker v1.0</span></div>",
        unsafe_allow_html=True)
    st.sidebar.divider()
    st.sidebar.markdown(
        "<div style='font-family:JetBrains Mono,monospace;font-size:10px;"
        "color:#6e7681;text-transform:uppercase;letter-spacing:0.1em;margin-bottom:8px'>FILTERS</div>",
        unsafe_allow_html=True)
    mutations_df = data["mutations"]
    date_range = None
    if not mutations_df.empty and "date" in mutations_df.columns:
        mutations_df["date"] = pd.to_datetime(mutations_df["date"])
        min_date = mutations_df["date"].min().date()
        max_date = mutations_df["date"].max().date()
        date_range = st.sidebar.date_input("date range", value=(min_date, max_date),
                                            min_value=min_date, max_value=max_date)
    region = "All"
    if not mutations_df.empty and "region" in mutations_df.columns:
        regions = ["All"] + sorted(mutations_df["region"].dropna().unique().tolist())
        region = st.sidebar.selectbox("region", regions)
    variant = "All"
    if not mutations_df.empty and "who_label" in mutations_df.columns:
        variants = ["All"] + sorted(mutations_df["who_label"].dropna().unique().tolist())
        variant = st.sidebar.selectbox("variant", variants)
    st.sidebar.divider()
    source_label = "sample data · 5k seqs"
    real_path = OUTPUT_DIR / "nextstrain_metadata.parquet"
    if real_path.exists():
        try:
            df_size = len(pd.read_parquet(real_path))
            if df_size > 5000:
                source_label = f"nextstrain · {df_size:,} seqs"
        except Exception:
            pass
    st.sidebar.markdown(
        f"<div style='font-family:JetBrains Mono,monospace;font-size:10px;"
        f"color:#3fb950;background:#0f2d1a;border:1px solid #1a4427;"
        f"border-radius:4px;padding:6px 10px'>● {source_label}</div>",
        unsafe_allow_html=True)
    return date_range, region, variant

def apply_filters(data, date_range, region, variant):
    df = data["mutations"].copy()
    if df.empty:
        return data
    df["date"] = pd.to_datetime(df["date"])
    if date_range and len(date_range) == 2:
        df = df[(df["date"] >= pd.Timestamp(date_range[0])) &
                (df["date"] <= pd.Timestamp(date_range[1]))]
    if region != "All" and "region" in df.columns:
        df = df[df["region"] == region]
    if variant != "All" and "who_label" in df.columns:
        df = df[df["who_label"] == variant]
    filtered = dict(data)
    filtered["mutations"] = df
    return filtered

def render_metrics(df):
    unique = set()
    for m in df["spike_mutations"].dropna():
        if isinstance(m, str):
            unique.update(json.loads(m))
    c1, c2, c3, c4, c5 = st.columns(5)
    with c1: st.metric("sequences", f"{len(df):,}")
    with c2: st.metric("unique mutations", f"{len(unique):,}")
    with c3: st.metric("who variants", df["who_label"].nunique() if "who_label" in df.columns else 0)
    with c4: st.metric("countries", df["country"].nunique() if "country" in df.columns else 0)
    with c5:
        avg = df["rbd_mut_count"].mean() if "rbd_mut_count" in df.columns else 0
        st.metric("avg rbd muts", f"{avg:.1f}")

def tab_overview(data, filtered):
    render_metrics(filtered["mutations"])
    st.divider()
    col_left, col_right = st.columns([3, 2])
    with col_left:
        section("variant dominance", "monthly frequency of WHO-labelled variants · 2020–2026")
        vf = data["variant_freq"]
        if not vf.empty:
            fig = dark_fig(timeline_plot(vf, mode="variant"))
            st.plotly_chart(fig, width='stretch', key="chart_1")
    with col_right:
        section("evolutionary events")
        events = [
            ("2020-02", "D614G",   "Transmissibility", "sweeps globally — universal backbone"),
            ("2020-09", "Alpha",   "Alpha",   "N501Y + P681H · 50% transmission advantage"),
            ("2020-10", "E484K",   "Immune",  "convergent immune escape · 3 independent lineages"),
            ("2021-04", "Delta",   "Delta",   "L452R + T478K · P681R furin cleavage"),
            ("2021-11", "Omicron", "Omicron", "30+ spike mutations · unprecedented immune evasion"),
            ("2022-10", "XBB",     "Omicron", "recombinant lineage · F486P ACE2 escape"),
            ("2023-08", "JN.1",    "Omicron", "BA.2.86 descendant · L455S ACE2 shift"),
        ]
        for date, label, var, desc in events:
            col = WHO_COLOURS.get(var, "#444")
            st.markdown(
                f"<div style='background:#161b22;border:1px solid #21262d;"
                f"border-left:3px solid {col};border-radius:0 6px 6px 0;"
                f"padding:7px 12px;margin:3px 0'>"
                f"<span style='font-family:JetBrains Mono,monospace;font-size:10px;color:#6e7681'>{date}</span><br>"
                f"<span style='font-weight:600;color:#e6edf3;font-size:13px'>{label}</span> "
                f"<span style='color:#8b949e;font-size:12px'>— {desc}</span></div>",
                unsafe_allow_html=True)

def tab_timeline(data):
    section("mutation frequency", "track how variants and mutations spread over time")
    mode = st.radio("view", ["variant dominance", "individual mutations"], horizontal=True)
    if mode == "individual mutations":
        all_muts = [m[0] for m in KEY_MUTATIONS]
        selected = st.multiselect("mutations", all_muts,
                                   default=["D614G", "N501Y", "E484K", "L452R", "P681H"])
        fig = dark_fig(timeline_plot(data["variant_freq"], mutation_freq=data["mutation_freq"],
                                      mode="mutation", selected_mutations=selected))
    else:
        fig = dark_fig(timeline_plot(data["variant_freq"], mode="variant"))
    st.plotly_chart(fig, width='stretch', key="chart_2")
    st.divider()
    section("selective sweep speed", "days from 5% → 50% global frequency")
    sweeps = data["sweeps"]
    if not sweeps.empty:
        col1, col2 = st.columns([2, 1])
        with col1:
            st.plotly_chart(dark_fig(sweep_speed_chart(sweeps)), key="chart_3", width='stretch')
        with col2:
            st.dataframe(
                sweeps[["mutation","first_seen","reached_50pct","sweep_days","peak_freq"]]
                      .rename(columns={"mutation":"mut","first_seen":"first seen",
                                       "reached_50pct":"50%","sweep_days":"days","peak_freq":"peak %"}),
                width='stretch', hide_index=True)

def tab_geographic(data):
    section("geographic spread", "variant origins and sequence contributions by country")
    geo = data["geo"]
    if geo.empty:
        st.info("no geographic data available")
        return
    st.plotly_chart(dark_fig(geographic_plotly(geo)), width='stretch', key="chart_4")
    col1, col2 = st.columns(2)
    with col1:
        section("first detection by variant")
        first = (geo.sort_values("first_detection").drop_duplicates(subset=["who_label"])
                    [["who_label","country","first_detection","sequences_count"]]
                    .rename(columns={"who_label":"variant","first_detection":"first detected","sequences_count":"seqs"}))
        st.dataframe(first, width='stretch', hide_index=True)
    with col2:
        section("top contributing countries")
        top = (geo.groupby("country")["sequences_count"].sum().reset_index()
                  .sort_values("sequences_count", ascending=False).head(12)
                  .rename(columns={"sequences_count":"sequences"}))
        st.dataframe(top, width='stretch', hide_index=True)

def tab_spike(data):
    section("spike protein", "1,273 amino acid map · mutational hotspots")
    hs = data["hotspots"]
    st.plotly_chart(dark_fig(spike_protein_map(hs)), width='stretch', key="chart_5")
    col1, col2 = st.columns(2)
    with col1:
        section("rbd zoom", "receptor binding domain · positions 306–527")
        st.plotly_chart(dark_fig(spike_rbd_detail(hs)), width='stretch', key="chart_6")
    with col2:
        section("mutation categories")
        import plotly.graph_objects as go
        from config import FUNCTIONAL_COLOURS
        cats = {}
        for _, _, cat, _ in KEY_MUTATIONS:
            cats[cat] = cats.get(cat, 0) + 1
        fig_pie = go.Figure(go.Pie(
            labels=list(cats.keys()), values=list(cats.values()),
            marker_colors=[FUNCTIONAL_COLOURS.get(c, "#888") for c in cats],
            hole=0.5, textinfo="label+percent",
            textfont=dict(size=10, color="#8b949e")))
        fig_pie.update_layout(showlegend=False, margin=dict(l=10,r=10,t=10,b=10),
                               paper_bgcolor="#0d1117", plot_bgcolor="#0d1117",
                               font=dict(color="#8b949e"))
        st.plotly_chart(fig_pie, width='stretch', key="chart_7")
    section("domain breakdown")
    if not hs.empty:
        ds = (hs.groupby("domain")
                .agg(mutations=("mutation_count","sum"), positions=("position","count"),
                     top=("position", lambda x: ", ".join(x.nlargest(3).astype(str))))
                .sort_values("mutations", ascending=False).reset_index())
        st.dataframe(ds, width='stretch', hide_index=True)

def tab_convergent(data):
    section("convergent evolution", "mutations that arose independently in multiple lineages")
    conv = data["convergent"]
    if conv.empty:
        st.info("no convergent evolution data available")
        return
    import plotly.express as px
    fig = px.bar(conv.nlargest(12, "appearance_count"), x="mutation", y="appearance_count",
                 color="appearance_count",
                 color_continuous_scale=[[0,"#1a4427"],[0.5,"#2ea043"],[1.0,"#3fb950"]],
                 labels={"mutation":"mutation","appearance_count":"lineage appearances"})
    fig.update_layout(coloraxis_showscale=False, showlegend=False,
                      xaxis_title="spike mutation", yaxis_title="independent appearances",
                      **PLOTLY_DARK)
    fig.update_traces(marker_line_width=0)
    st.plotly_chart(fig, width='stretch', key="chart_8")
    st.markdown(
        "<div style='background:#161b22;border:1px solid #21262d;border-left:3px solid #3fb950;"
        "border-radius:0 6px 6px 0;padding:10px 14px;font-size:12px;color:#8b949e;margin:8px 0'>"
        "When the same mutation appears independently across unrelated lineages, it signals strong "
        "selective pressure. <strong style='color:#e6edf3'>E484K</strong> is the textbook case — "
        "it evolved separately in Beta (South Africa), Gamma (Brazil), and Iota (USA) with no shared ancestor."
        "</div>", unsafe_allow_html=True)
    col1, col2 = st.columns(2)
    with col1:
        section("top convergent mutations")
        st.dataframe(
            conv.nlargest(8,"appearance_count")[["mutation","appearance_count","lineages"]]
                .rename(columns={"mutation":"mutation","appearance_count":"variants","lineages":"found in"}),
            width='stretch', hide_index=True)
    with col2:
        section("what this means")
        st.markdown(
            "<div style='font-size:12px;color:#8b949e;line-height:1.8'>"
            "Each bar = a mutation the virus 'discovered' multiple times independently.<br><br>"
            "High bars = strong selection pressure.<br><br>"
            "These are almost always linked to <span style='color:#58a6ff'>ACE2 binding</span> "
            "or <span style='color:#f85149'>immune escape</span>.</div>",
            unsafe_allow_html=True)

def tab_data(filtered):
    section("sequence database", "raw data · searchable · exportable")
    df = filtered["mutations"].copy()
    if df.empty:
        st.info("no sequences match current filters")
        return
    search = st.text_input("", placeholder="search by sample ID, country, or mutation")
    if search:
        mask = (df["strain"].str.contains(search, case=False, na=False) |
                df["country"].str.contains(search, case=False, na=False) |
                df["spike_mutations"].str.contains(search, case=False, na=False))
        df = df[mask]
    show_cols = [c for c in ["strain","date","country","region","lineage",
                              "who_label","spike_mut_count","rbd_mut_count"] if c in df.columns]
    st.dataframe(
        df[show_cols].rename(columns={
            "strain":"sample id","date":"date","country":"country","region":"region",
            "lineage":"pango","who_label":"who","spike_mut_count":"spike muts","rbd_mut_count":"rbd muts"})
        .sort_values("date", ascending=False).head(500),
        width='stretch', hide_index=True)
    st.markdown(
        f"<div style='font-family:JetBrains Mono,monospace;font-size:10px;color:#6e7681;"
        f"margin-top:6px'>showing {min(500,len(df)):,} of {len(df):,} sequences</div>",
        unsafe_allow_html=True)
    st.download_button("↓ export csv", df[show_cols].to_csv(index=False),
                        "sars_cov2_mutations.csv", "text/csv")

def tab_about():
    section("about")
    st.markdown("""
<div style='font-size:13px;color:#8b949e;line-height:1.8;max-width:700px'>
<p>A genomic surveillance pipeline tracking how SARS-CoV-2 spike protein mutations emerge,
spread geographically, and affect viral fitness. Uses real sequence data from the Nextstrain
open database — the same source used by the WHO and public health agencies worldwide.</p>
<br>
<table style='border-collapse:collapse;width:100%;font-family:JetBrains Mono,monospace;font-size:11px'>
<tr style='border-bottom:1px solid #21262d'><td style='padding:8px 16px 8px 0;color:#6e7681'>fetch_data.py</td><td style='padding:8px 0;color:#8b949e'>downloads MN908947 reference · fetches Nextstrain metadata</td></tr>
<tr style='border-bottom:1px solid #21262d'><td style='padding:8px 16px 8px 0;color:#6e7681'>alignment.py</td><td style='padding:8px 0;color:#8b949e'>parses amino-acid substitutions · builds mutation table</td></tr>
<tr style='border-bottom:1px solid #21262d'><td style='padding:8px 16px 8px 0;color:#6e7681'>frequency_analysis.py</td><td style='padding:8px 0;color:#8b949e'>variant frequency · sweep speed · convergent evolution · hotspots</td></tr>
<tr style='border-bottom:1px solid #21262d'><td style='padding:8px 16px 8px 0;color:#6e7681'>plots.py</td><td style='padding:8px 0;color:#8b949e'>all plotly figures · choropleth · spike protein map</td></tr>
<tr><td style='padding:8px 16px 8px 0;color:#6e7681'>app.py</td><td style='padding:8px 0;color:#8b949e'>this streamlit dashboard</td></tr>
</table>
<br>
<p style='color:#6e7681;font-size:11px;font-family:JetBrains Mono,monospace'>
data: nextstrain.org · gisaid.org · ncbi.nlm.nih.gov/nuccore/MN908947<br>
sequences gratefully acknowledged from all contributing laboratories worldwide.
</p></div>""", unsafe_allow_html=True)

def main():
    st.markdown(
        "<div style='padding:8px 0 4px'>"
        "<span style='font-family:JetBrains Mono,monospace;font-size:22px;"
        "font-weight:600;color:#e6edf3;letter-spacing:-0.02em'>SARS-CoV-2</span> "
        "<span style='font-family:JetBrains Mono,monospace;font-size:22px;"
        "font-weight:300;color:#58a6ff'>// mutation tracker</span><br>"
        "<span style='font-family:JetBrains Mono,monospace;font-size:11px;"
        "color:#6e7681'>spike protein evolution · 2020–2026 · nextstrain open data</span>"
        "</div>", unsafe_allow_html=True)
    st.divider()
    with st.spinner("loading..."):
        data = load_all_data()
    date_range, region, variant = render_sidebar(data)
    filtered = apply_filters(data, date_range, region, variant)
    tabs = st.tabs(["overview","mutation timeline","geographic spread",
                    "spike protein","convergent evolution","data explorer","about"])
    with tabs[0]: tab_overview(data, filtered)
    with tabs[1]: tab_timeline(data)
    with tabs[2]: tab_geographic(data)
    with tabs[3]: tab_spike(data)
    with tabs[4]: tab_convergent(data)
    with tabs[5]: tab_data(filtered)
    with tabs[6]: tab_about()

if __name__ == "__main__":
    main()
